#include "uncolfalgo.h"

#include "ChiPhysics/FieldFunction/fieldfunction.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Executes the algorithm.*/
//void UncolFAlgo::Solver::Execute()
//{
//  //======================================== Determine rays and angles
//  int SV_counter=0;
//  for (const auto& SV : SVs)
//  {
//    double total_solid_angle = 0.0;
//    chi_log.Log(LOG_0) << "SV " << ++SV_counter << " of " << SVs.size();
//
//    for (const auto& IS : ISs)
//    {
//      auto IS_point = GetCentroidFromList(IS);
//
//      auto   ray_direction   = (IS_point - SV.SP).Normalized();
//      double ray_solid_angle = ComputeSolidAngle(SV.SP, IS, geometry_type);
//
//      Ray ray(SV,                       //Source Volume
//              ray_direction,            //Direction
//              ray_solid_angle,          //Solid angle
//              SV.owning_cell_global_id, //Owning cell global id
//              SV.SP);                   //Initial position
//
//
//      bool alive = true; while (alive) alive = RayTrace(ray);
//
//      total_solid_angle += ray_solid_angle;
//    }//for IS_face
//
//    chi_log.Log(LOG_0)
//      << "SV " << SV_counter << " of " << SVs.size() << " "
//      << "Total solid angle: " << total_solid_angle;
//  }//for SP
//
//  ComputePWLDTransformations();
//
//  std::vector<chi_physics::FieldFunction*> ff_list;
//  for (int m=0; m<num_moments; ++m)
//  {
//    auto ff = new chi_physics::FieldFunction(
//      "asd",
//      chi_physics::FieldFunctionType::DFEM_PWL,
//      pwl,
//      &uk_man,
//      0,
//      &phi_new_local);
//
//    ff_list.push_back(ff);
//  }
//
//  chi_physics::FieldFunction::ExportMultipleFFToVTK("UPhi",ff_list);
//
//}

//###################################################################
/**Executes the algorithm.*/
void UncolFAlgo::Solver::Execute()
{
  //======================================== Build schedule of stages
  int num_SVs = SVs.size();
  int num_IS_levels = IS_levels.size();

  typedef std::vector<int> VecInt;
  typedef std::vector<VecInt> MatInt;

  VecInt SV_level_inds(num_SVs,-1);
  int num_stages=0;
  {
    int cur_index=-1;
    for (auto& SV_level_ind : SV_level_inds)
      SV_level_ind = cur_index--;

    num_stages = std::abs(cur_index) - 2 + num_IS_levels;
    chi_log.Log(LOG_0) << "Number of stages to be executed: " << num_stages;
  }

  //======================================== Execute stages
  for (int stage=0; stage<num_stages; ++stage)
  {
    chi_log.Log(LOG_0) << "Stage " << stage;
    for (int i=0; i<num_SVs; ++i)
    {
      int& ell = SV_level_inds[i] += 1;
      auto& SV = SVs[i];

      if (not (ell >=0 and ell < num_IS_levels)) continue;

      chi_log.Log(LOG_0) << "  SV " << i << " level " << ell << std::endl;

      int num_ISpoints = IS_levels[ell].cone_directions.size();
      for (int q=0; q<num_ISpoints; ++q)
      {
        auto   ray_direction   = IS_levels[ell].cone_directions[q];
        double ray_solid_angle = IS_levels[ell].cone_solid_angles[q];

        Ray ray(ray_direction,            //Direction
                ray_solid_angle,          //Solid angle
                SV.owning_cell_global_id, //Owning cell global id
                SV.SP,                    //Initial position
                ell,                      //Reference IS_level
                q,                        //Reference ISP_index
                i);                       //Source Volume index

        bool alive = true; while (alive) alive = RayTrace(ray);
      }//for q
    }
  }

  ComputePWLDTransformations();

  std::vector<chi_physics::FieldFunction*> ff_list;
  for (int m=0; m<num_moments; ++m)
  {
    auto ff = new chi_physics::FieldFunction(
      "asd",
      chi_physics::FieldFunctionType::DFEM_PWL,
      pwl,
      &uk_man,
      0,
      &phi_new_local);

    ff_list.push_back(ff);
  }

  chi_physics::FieldFunction::ExportMultipleFFToVTK("UPhi",ff_list);

}