#include "uncolfalgo.h"

#include "ChiMesh/Cell/cell_polyhedron.h"
#include "ChiMath/SpatialDiscretization/PiecewiseLinear/CellViews/pwl_polyhedron.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Process a track inside of the designated source volume.*/
void UncolFAlgo::Solver::ProcessTrackInsideSV(Ray& ray,
                                              const chi_mesh::Vector3& ray_destination,
                                              chi_mesh::Cell& cell,
                                              double sig_t)
{
  auto& SV = SVs[ray.SV_index];

  //======================================== Grab cell fe info and
  //                                         quadrature information
  auto cell_pwl_view = pwl->MapFeViewL(cell.local_id);
  const int NDOFS = cell_pwl_view->dofs;

  if (geometry_type == GeometryType::THREED_CARTESIAN)
  {
    auto polyh_cell = (chi_mesh::CellPolyhedron*)&cell;
    PolyhedronFEView polyh_pwl_view(polyh_cell,grid,pwl);
    polyh_pwl_view.PreCompute();

    //====================================== Compute cell volume
    double V = 0.0;
    {
      auto& v0 = *grid->vertices[SV.ref_cell->vertex_ids[0]];
      auto& v1 = *grid->vertices[SV.ref_cell->vertex_ids[1]];
      auto& v2 = *grid->vertices[SV.ref_cell->vertex_ids[2]];
      auto& v3 = *grid->vertices[SV.ref_cell->vertex_ids[3]];

      //====================================== Build jacobian
      auto v01 = v1-v0;
      auto v02 = v2-v0;
      auto v03 = v3-v0;

      chi_mesh::Matrix3x3 J;
      J.SetColJVec(0,v01);
      J.SetColJVec(1,v02);
      J.SetColJVec(2,v03);


      V = J.Det()/6.0;
    }//cell V

    //====================================== Compute projected area
    double A_p = 0.0;
    for (const auto& face : SV.faces)
    {
      auto v01 = face[1] - face[0];
      auto v02 = face[2] - face[0];

      double area = 0.5*v01.Cross(v02).Norm();
      auto normal = v01.Cross(v02).Normalized();

      if (normal.Dot(ray.omega)>0.0)
        A_p += area*normal.Dot(ray.omega);
    }

    //====================================== Compute tet averaged psi
    double Lmax = 3.0*V/A_p;
    double tau_max = sig_t*Lmax;

    double psi  = 6.0*tau_max - 3.0*tau_max*tau_max + tau_max*tau_max*tau_max -
                  6.0 + 6.0*exp(-tau_max);
           psi /= tau_max*tau_max*tau_max;
           psi *= SV.magnitude/sig_t;

    double Ylm = 1.0;
    double w_psi_Ylm = ray.solid_angle*psi*Ylm;

    std::vector<double> IntOmega_IntV_PsiShape_i(NDOFS,0.0);

    for (auto& face_data : polyh_pwl_view.face_data)
      for (auto& side_data : face_data.sides)
        for (int i=0; i<NDOFS; ++i)
          for (double shape_i_qp : side_data.qp_data[i].shape_qp)
            IntOmega_IntV_PsiShape_i[i] +=
              w_psi_Ylm*shape_i_qp*side_data.detJ*(1.0/6.0);


    for (int i=0; i<NDOFS; ++i)
    {
      int ir = pwl->MapDFEMDOFLocal(&cell,i,&uk_man,0,0);
      phi_old_local_src[ir] += IntOmega_IntV_PsiShape_i[i];
    }
  }//3d cartesian
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "UncolFAlgo::Solver::ProcessTrackInsideSV."
      << "Unsupported geometry type.";
    exit(EXIT_FAILURE);
  }

}