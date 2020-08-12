#include "uncolfalgo.h"

#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/Cell/cell_polygon.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/** Initializes.*/
void UncolFAlgo::Solver::Initialize()
{
  //============================================= Get reference to grid
  auto cur_handler = chi_mesh::GetCurrentHandler();
  grid = cur_handler->GetGrid();

  //============================================= Determine geometry type
  if (grid->local_cells[0].Type() == chi_mesh::CellType::SLAB)
    geometry_type = GeometryType::ONED_SLAB;
  if (grid->local_cells[0].Type() == chi_mesh::CellType::POLYGON)
    geometry_type = GeometryType::TWOD_CARTESIAN;
  if (grid->local_cells[0].Type() == chi_mesh::CellType::POLYHEDRON)
    geometry_type = GeometryType::THREED_CARTESIAN;

  //============================================= Add spatial discretization
  chi_log.Log(LOG_0) << "Computing cell matrices";
  pwl = new SpatialDiscretization_PWL(2);
  pwl->AddViewOfLocalContinuum(grid);

  //============================================= Add unknown management
  for (int m=0; m<num_moments; ++m)
  {
    uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N,num_groups);
    for (int g=0; g<num_groups; ++g)
      uk_man.SetUnknownComponentTextName(m,g,
                                         "m"+std::to_string(m)+
                                         "_g"+std::to_string(g));
  }


  //============================================= Determine dof ordering
  pwl->OrderNodesDFEM(grid);
  local_dof_count = pwl->GetNumLocalDOFs(grid,&uk_man);
  globl_dof_count = pwl->GetNumLocalDOFs(grid,&uk_man);

  phi_old_local.resize(local_dof_count,0.0);
  phi_old_local_src.resize(local_dof_count,0.0);
  phi_new_local.resize(local_dof_count,0.0);

  //============================================= Initialize cell Ak
  Ak.resize(grid->local_cells.size());
  for (auto& cell : grid->local_cells)
  {
    int Ndofs = cell.vertex_ids.size();
    Ak[cell.local_id] = MatDbl(Ndofs,VecDbl(Ndofs,0.0));
  }

  int IS_refinement_level = 1;
  int SV_refinement_level = 0;

  //============================================= Build Interpolation Surfaces (ISs)
  chi_log.Log(LOG_0) << "Building Interpolation Surfaces";
  for (auto& cell : grid->local_cells)
    for (auto& face : cell.faces)
      if ((not face.IsNeighborLocal(grid)) or face.neighbor < 0)
      {
        if (geometry_type == GeometryType::ONED_SLAB)
        {
          VertList vertices(1);
          vertices[0] = face.centroid;

          ISs.push_back(vertices);
        }//1D
        else if (geometry_type == GeometryType::TWOD_CARTESIAN)
        {
          auto& v0 = *grid->vertices[face.vertex_ids[0]];
          auto& v1 = *grid->vertices[face.vertex_ids[1]];
          auto v01 = (v1-v0);

          double dz = 1.0/(IS_refinement_level + 1);
          for (int l=0; l < (IS_refinement_level + 1); ++l)
          {
            double d0 = l*dz;
            double d1 = (l+1)*dz;

            VertList vertices(2);
            vertices[0] = v0 + d0*v01;
            vertices[1] = v0 + d1*v01;

            ISs.push_back(vertices);
          }

        }//2D
        else if (geometry_type == GeometryType::THREED_CARTESIAN)
        {
          const int vlast = face.vertex_ids.size()-1;
          auto& v0 = face.centroid;
          for (int v=0; v<face.vertex_ids.size(); ++v)
          {
            auto& v1 = *grid->vertices[face.vertex_ids[v]];
            auto& v2 = *grid->vertices[face.vertex_ids[(v<vlast)? v+1 : 0]];

            auto triangles = SubdivideTriangle({{v0,v1,v2}},IS_refinement_level);

            for (const auto& triangle : triangles)
              ISs.push_back(triangle);
          }
        }//3D
      }//if face is interface

  //============================================= Build Source Volumes (SVs)
  chi_log.Log(LOG_0) << "Building Source Volumes";
  for (auto& cell : grid->local_cells)
  {
    if (cell.material_id == 1) //TODO: Rewordk, if material has source
    {
      double sig_t = 1.0; //TODO: Properly get this from material

      if (geometry_type == GeometryType::ONED_SLAB)
      {
        auto& v0 = *grid->vertices[cell.vertex_ids[0]];
        auto& v1 = *grid->vertices[cell.vertex_ids[1]];

        VertCollection faces(2);
        faces[0] = {v0};
        faces[1] = {v1};

        SVs.emplace_back(cell.centroid, 1.0, cell.global_id, faces, sig_t);
      }
      else if (geometry_type == GeometryType::TWOD_CARTESIAN)
      {
        auto& v0 = cell.centroid;
        grid->vertices.push_back(new chi_mesh::Vertex(v0));
        for (auto& face : cell.faces)
        {
          auto& v1 = *grid->vertices[face.vertex_ids[0]];
          auto& v2 = *grid->vertices[face.vertex_ids[1]];
          
          auto triangles = SubdivideTriangle({{v0,v1,v2}}, SV_refinement_level);
          
          for (const auto& triangle_verts : triangles)
          {
            VertCollection faces = {{triangle_verts[0], triangle_verts[1]},
                                    {triangle_verts[1], triangle_verts[2]},
                                    {triangle_verts[2], triangle_verts[0]}};

            auto vc = GetCentroidFromList(triangle_verts);

            auto new_cell = SpawnTriangle(triangle_verts,grid);

            SVs.emplace_back(vc, 1.0, cell.global_id,faces, sig_t);

            SVs.back().ref_cell=new_cell;
          }//for triangle
        }//for face
      }//polygon
      else if (geometry_type == GeometryType::THREED_CARTESIAN)
      {
        for (auto& face : cell.faces)
        {
          auto& v0 = cell.centroid;
          for (int v=0; v<face.vertex_ids.size(); v++)
          {
            int V = face.vertex_ids.size()-1;

            auto& v1 = face.centroid;
            auto& v2 = *grid->vertices[face.vertex_ids[v]];
            auto& v3 = *grid->vertices[face.vertex_ids[(v<V)?v+1:0]];

            auto tets = SubDivideTetrahedron({{v0,v1,v2,v3}},SV_refinement_level);

            for (const auto& tet : tets)
            {
              VertCollection faces = {{tet[0],tet[2],tet[1]},
                                      {tet[2],tet[3],tet[1]},
                                      {tet[3],tet[0],tet[1]},
                                      {tet[0],tet[3],tet[2]}};

              auto SP = GetCentroidFromList(tet);

              SourceVolume source_volume(SP, 1.0, cell.global_id, faces, sig_t);
              source_volume.ref_cell = SpawnTetrahedron(tet,grid);

              SVs.push_back(source_volume);
            }//for tet
          }//for v
        }//for face
      }//polygon
    }//material has source
  }



}