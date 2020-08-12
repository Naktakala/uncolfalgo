#include "uncolfalgo.h"

#include "ChiMesh/Raytrace/raytracing.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Traces a ray over cells.*/
bool UncolFAlgo::Solver::RayTrace(Ray &ray)
{
  auto& cell = *grid->cells[ray.current_cell];

  //======================================== Process material
  double sig_t = 1.0;

  //======================================== Raytrace to surface
  chi_mesh::Vector3 surface_intersection_point;
  double d_to_surface = 1.0e3;

  chi_mesh::RayDestinationInfo ray_dest_info =
    chi_mesh::RayTrace(
      grid,                         //[Input] Grid
      &cell,                        //[Input] Current cell
      ray.current_position,         //[Input] Current position
      ray.omega,                    //[Input] Current direction
      d_to_surface,                 //[Otput] Distance to next surface
      surface_intersection_point    //[Otput] Intersection point at next surf
      );

  //======================================== Process Track
  if (ray.current_cell == ray.SV.owning_cell_global_id)
  {
    ProcessTrackInsideSV(ray,surface_intersection_point, cell, sig_t);
    ray.exit_point = surface_intersection_point;
  }
  else
  {
    ProcessTrackOutsideSV(ray, surface_intersection_point, cell, sig_t);
    ray.tau  += d_to_surface*sig_t;
  }

  //======================================== Advance the ray
  ray.current_position = surface_intersection_point;
  ray.current_cell     = ray_dest_info.destination_face_neighbor;

  //======================================== Determine if ray is still alive
  bool ray_alive = true;
  if (ray_dest_info.destination_face_neighbor < 0) ray_alive = false;

  return ray_alive;
}