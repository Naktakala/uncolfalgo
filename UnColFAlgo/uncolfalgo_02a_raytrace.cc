#include "uncolfalgo.h"

#include "ChiMesh/Raytrace/raytracing.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Traces a ray over cells.*/
bool UncolFAlgo::Solver::RayTrace(Ray &ray)
{
//  chi_log.Log() << "Raytracer";
  auto& cell = *grid->cells[ray.current_cell];
  LDFEIS& IS = IS_levels[ray.IS_level];
  auto& IS_tau = IS.SPig_taus[ray.SV_index];
  auto& SV = SVs[ray.SV_index];

  //======================================== Process material
  double sig_t = 1.0;

  //======================================== Raytrace to surface
  chi_mesh::Vector3 surface_intersection_point;
  double d_to_surface = 1.0e3;

  chi_mesh::RayDestinationInfo ray_dest_info =
    chi_mesh::RayTrace(
      grid,                         //[Input]  Grid
      &cell,                        //[Input]  Current cell
      ray.current_position,         //[Input]  Current position
      ray.omega,                    //[Input]  Current direction
      d_to_surface,                 //[Output] Distance to next surface
      surface_intersection_point    //[Output] Intersection point at next surf
      );

  //========================================
  double R_level_i_sqr = IS.radius_inner_sqr;
  double R_level_o_sqr = IS.radius_outer_sqr;
  double R_i_sqr = (ray.current_position - SV.SP).NormSquare();
  double R_o_sqr = (surface_intersection_point - SV.SP).NormSquare();

  //======================================== Process Track inside SV
  if (ray.current_cell == SV.owning_cell_global_id)
  {
//    if (R_i_sqr >= R_level_i_sqr and R_i_sqr <= R_level_o_sqr)
//      ProcessTrackInsideSV(ray,surface_intersection_point, cell, sig_t);
    ray.exit_point = surface_intersection_point;
  }
  //======================================== Process Track outside SV
  else
  {
    auto track_endpoint = surface_intersection_point;
    double d_to_endpoint = d_to_surface;

    //================================= If within the reference IS_level
    if (R_i_sqr >= R_level_i_sqr and R_i_sqr <= R_level_o_sqr)
    {
      //===================== If intersecting the IS
      // This means it needs to terminate
      // exactly on the IS and set the ISP's
      // reference tau value.
      if (R_o_sqr > R_level_o_sqr)
      {
        double d_to_ISo = sqrt(R_level_o_sqr) - sqrt(R_i_sqr);
        d_to_endpoint = d_to_ISo;
        track_endpoint = ray.current_position + d_to_ISo * ray.omega;
        IS_tau[ray.ISP_index] = ray.tau + d_to_ISo * sig_t;
      }
    }

    //================================= If ahead of a child IS_level
    if (IS.parent != nullptr and
        R_i_sqr < R_level_i_sqr and R_o_sqr >= R_level_i_sqr)
    {
      double d_to_ISi = sqrt(R_level_i_sqr) - sqrt(R_i_sqr);
      ray.current_position += d_to_ISi * ray.omega;

      int ref_SQ_index = IS.parent_SSQ_mapping[ray.ISP_index];

      ray.tau = IS.parent->SPig_taus[ray.SV_index][ref_SQ_index];
      d_to_endpoint = d_to_surface - d_to_ISi;
    }

//    ProcessTrackOutsideSV(ray, track_endpoint, cell, sig_t);
    ray.tau  += d_to_endpoint*sig_t;
  }

  //======================================== Advance the ray
  ray.current_position = surface_intersection_point;
  ray.current_cell     = ray_dest_info.destination_face_neighbor;

  //======================================== Determine if ray is still alive
  bool ray_alive = true;
  if (ray_dest_info.destination_face_neighbor < 0) ray_alive = false;

  //======================================== Limit extent
  if (R_o_sqr > R_level_o_sqr) ray_alive = false;

//  chi_log.Log() << "exit Raytracer";
  return ray_alive;
}