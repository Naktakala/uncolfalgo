#include "uncolfalgo.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Executes the algorithm.*/
void UncolFAlgo::Solver::Execute()
{
  std::vector<Ray> rays;

  //======================================== Determine rays and angles
  int SV_counter=0;
  for (const auto& SV : SVs)
  {
    double total_solid_angle = 0.0;
    ++SV_counter;
    chi_log.Log(LOG_0)
      << "SV " << SV_counter << " of " << SVs.size();
    for (const auto& IS : ISs)
    {
      auto IS_point = GetCentroidFromList(IS);

      chi_mesh::Vector3 ray_direction = (IS_point - SV.SP).Normalized();
      double          ray_solid_angle = ComputeSolidAngle(SV.SP, IS, geometry_type);

      Ray ray(SV,                       //Source Volume
              ray_direction,            //Direction
              ray_solid_angle,          //Solid angle
              SV.owning_cell_global_id, //Owning cell global id
              SV.SP);                   //Initial position


      bool alive = true; while (alive) alive = RayTrace(ray);

      total_solid_angle += ray.solid_angle;
    }//for IS_face

//    chi_log.Log(LOG_0)
//      << "SV " << SV_counter << " of " << SVs.size() << " "
//      << "Total solid angle: " << total_solid_angle;
  }//for SP

  ComputePWLDTransformations();

  for (int m=0; m<num_moments; ++m)
  {

  }
  auto ff = new chi_physics::FieldFunction(
    "asd",
    chi_physics::FieldFunctionType::DFEM_PWL,
    pwl,
    &uk_man,
    0,
    &phi_new_local);

  std::vector<chi_physics::FieldFunction*> ff_list = {ff};

  chi_physics::FieldFunction::ExportMultipleFFToVTK("UPhi",ff_list);

}