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

  //============================================= Initialize cell Ak-matrices
  Ak.resize(grid->local_cells.size());
  for (auto& cell : grid->local_cells)
  {
    int Ndofs = cell.vertex_ids.size();
    Ak[cell.local_id] = MatDbl(Ndofs,VecDbl(Ndofs,0.0));
  }

  BuildSVs();
  BuildISs();

}