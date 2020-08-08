#include "uncolfalgo.h"

#include "ChiMesh/Raytrace/raytracing.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Traces a ray over cells.*/
bool UncolFAlgo::Solver::RayTrace(Ray &ray)
{
  auto& cell = *grid->cells[ray.current_cell];
//  std::cout << "Cell " << cell.global_id << " " << ray.omega.PrintS() << "\n";

  //======================================== Process material
  double sig_t = 1.0;

  //======================================== Raytrace to surface
  chi_mesh::Vector3 surface_intersection_point;
  double d_to_surface = 1.0e15;

  chi_mesh::RayDestinationInfo ray_dest_info =
    chi_mesh::RayTrace(
      grid,                         //[Input] Grid
      &cell,                        //[Input] Current cell
      ray.current_position,         //[Input] Current position
      ray.omega,                    //[Input] Current direction
      d_to_surface,                 //[Otput] Distance to next surface
      surface_intersection_point);  //[Otput] Intersection point at next surf

  //======================================== Contribute to Y
  ContributeToY(ray,surface_intersection_point,cell,sig_t);

  //======================================== Grab exit point E if applicable
  if (ray.current_cell == ray.SV.owning_cell_global_id)
    ray.exit_point = surface_intersection_point;
  else
    ray.tau  += d_to_surface*sig_t;

  //======================================== Advance the ray
  ray.current_position = surface_intersection_point;
  ray.current_cell     = ray_dest_info.destination_face_neighbor;

  //======================================== Determine if ray is still alive
  bool ray_alive = true;
  if (ray_dest_info.destination_face_neighbor < 0) ray_alive = false;

  return ray_alive;
}