#include "uncolfalgo.h"

#include "ChiMesh/Raytrace/raytracing.h"

//###################################################################
/**Process a track outside of the designated source volume.*/
void UncolFAlgo::Solver::ProcessTrackOutsideSV(Ray& ray,
                                               const chi_mesh::Vector3& ray_destination,
                                               chi_mesh::Cell& cell,
                                               double sig_t)
{
  if (cell.global_id == ray.SV.owning_cell_global_id) return;

  //======================================== Functor for raytracing
  struct RAY_TRACER
  {
    chi_mesh::MeshContinuum* grid;
    chi_mesh::Cell& cell;

    RAY_TRACER(chi_mesh::MeshContinuum* in_grid,
               chi_mesh::Cell& in_cell) :
      grid(in_grid), cell(in_cell)
    {}

    chi_mesh::Vector3 operator()(const chi_mesh::Vector3& from_point,
                    const chi_mesh::Vector3& direction)
    {
      chi_mesh::Vector3 intersection_point;
      double d_to_surface = 1.0e5;

      chi_mesh::RayDestinationInfo ray_dest_info =
        chi_mesh::RayTrace(
          grid,                         //[Input] Grid
          &cell,                        //[Input] Current cell
          from_point,                   //[Input] Current position
          direction,                    //[Input] Current direction
          d_to_surface,                 //[Otput] Distance to next surface
          intersection_point,           //[Otput] Intersection point at next surf
          1.0e-8,
          1.0e-10
          );


      return intersection_point;
    }
  };

  RAY_TRACER RayTracerCurCell(grid, cell);
  RAY_TRACER RayTracerSrcCell(grid, *ray.SV.ref_cell);
  double sig_t_s = ray.SV.sig_t;

  //======================================== Overall track length
  double tracklength = (ray_destination - ray.current_position).Norm();

  //======================================== Populate segment lengths
  segment_lengths.clear();
  segment_lengths.push_back(tracklength);

  chi_mesh::PopulateRaySegmentLengths(grid,&cell,
                                      segment_lengths,
                                      ray.current_position,
                                      ray_destination,
                                      ray.omega);

  //======================================== Grab cell fe info and
  //                                         quadrature information
  auto cell_pwl_view = pwl->MapFeViewL(cell.local_id);
  const int NDOFS = cell_pwl_view->dofs;

  auto& qpoints = legendre_quadrature.abscissae;
  auto& qweights = legendre_quadrature.weights;

  //======================================== Process each segment
  auto& E = ray.exit_point;
  auto& I = ray.current_position;
  double s_i = (I-E).Norm();
  for (auto segment_length : segment_lengths)
  {
    double s_f = s_i + segment_length;

    //============================ Segment start & end position
    auto r_i = E + s_i*ray.omega;
    auto r_f = E + s_f*ray.omega;

    //============================ Shape values at segment start & end
    cell_pwl_view->ShapeValues(r_i,shape_values_i);
    cell_pwl_view->ShapeValues(r_f,shape_values_f);

    int q_last = qpoints.size()-1;
    for (int q=0; q<qpoints.size(); ++q)
    {
      double xq = (qpoints[q]+1.0)/2.0;                          //x at q
      double xq_m1 = (q==0     )? 0.0 : (qpoints[q-1]+1.0)/2.0; //x at q-1
      double xq_p1 = (q==q_last)? 0.0 : (qpoints[q+1]+1.0)/2.0; //x at q+1

      double xq_mh = (q==0     )? 0.0 : 0.5*(xq+xq_m1);         //x at q-half
      double xq_ph = (q==q_last)? 1.0 : 0.5*(xq+xq_p1);         //x at q+half

      double sq    = s_i + xq*segment_length; //s at q
      double sq_mh = s_i + xq_mh*segment_length; //s at q-half
      double sq_ph = s_i + xq_ph*segment_length; //s at q+half

      double beta_q  = exp(-1.0*sig_t*(sq_mh - sq)) -
                       exp(-1.0*sig_t*(sq_ph - sq));
             beta_q /= sig_t;

      auto  Q = E + sq*ray.omega;

      double EI_norm = (I-E).Norm();
      double sig_t_outside = (std::fabs(EI_norm) > 1.0e-10)? ray.tau/EI_norm : 0.0;

      //=============================== Projection on source volume faces
      double J_q = 0.0;
      for (int f=0; f<ray.SV.faces.size(); ++f)
      {
        const auto& face = ray.SV.ref_cell->faces[f];
        const auto& face_verts = ray.SV.faces[f];

        if (face.normal.Dot(face_verts.front()-Q)<=0.0) continue;

        double domega_SC = ComputeSolidAngle(Q,face_verts,geometry_type);

        auto CC = GetCentroidFromList(face_verts);
        auto omega_QC = (CC - Q).Normalized();

        auto I_prime = RayTracerCurCell(Q, omega_QC);
        auto C       = RayTracerSrcCell(CC, -1.0 * omega_QC);

        double tau_cq = sig_t_outside*(I_prime - C).Norm() +
                        sig_t*(I_prime - Q).Norm();

        double L = (ray.SV.SP - E).Norm();

        double psi_exit = ray.SV.magnitude*(1.0-exp(-sig_t_s*L))/sig_t_s;

        double psi_q = psi_exit*exp(-tau_cq);

        double Ylm = 1.0;

        J_q += domega_SC * Ylm * psi_q;
      }//for subcone

      for (int i=0; i<NDOFS; ++i)
      {
        double bi = shape_values_i[i]*(1.0-xq) + shape_values_f[i]*(xq);

        for (int j=0; j<NDOFS; ++j)
        {
          double bj = shape_values_i[j]*(1.0-xq) + shape_values_f[j]*(xq);

          Ak[cell.local_id][i][j] +=
            ray.solid_angle*0.5*segment_length*qweights[q]*sq*sq*bi*bj;
        }//for j

        double Iq = beta_q*sq*sq*J_q;

        int ir = pwl->MapDFEMDOFLocal(&cell,i,&uk_man,0,0);
        phi_old_local[ir] += ray.solid_angle*bi*Iq;
      }//for i
    }//for q

    s_i = s_f;
  }
}