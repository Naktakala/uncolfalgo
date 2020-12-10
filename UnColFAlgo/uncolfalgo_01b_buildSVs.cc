#include "uncolfalgo.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/**Builds the Source Volumes.*/
void UncolFAlgo::Solver::BuildSVs()
{
  int SV_refinement_level = 0;

  //============================================= Build local Source Volumes
  chi_log.Log(LOG_0) << "Building Source Volumes";
  for (auto& cell : grid->local_cells)
  {
    if (cell.material_id == 1) //TODO: Rework, if material has source
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
  }//cell

  //============================================= Consolidating SV's
  int num_local_SVs = SVs.size();
  std::vector<int> locI_num_local_SVs(chi_mpi.process_count,0);

  MPI_Allgather(&num_local_SVs,             //sendbuf
                1,                          //sendcnt
                MPI_INT,                    //sendtype
                locI_num_local_SVs.data(),  //recvbuf
                1,                          //recvcnt
                MPI_INT,                    //recvtype
                MPI_COMM_WORLD);            //comm


}