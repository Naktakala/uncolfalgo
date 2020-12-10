#include "uncolfalgo.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Builds the Interpolation Surfaces.*/
void UncolFAlgo::Solver::BuildISs()
{

  int IS_refinement_level = 4;

  chi_log.Log(LOG_0) << "Building Interpolation Surfaces";


  LDFEISs_temp.GenerateInitialRefinement(IS_refinement_level);

  std::vector<LDFELocalRefinement> NO_LOCREF;

  int num_SVs_G = SVs.size()*1;
  IS_levels.reserve(3);
  IS_levels.emplace_back(3, NO_LOCREF, 3.0, -1.0, num_SVs_G, nullptr);
  IS_levels.emplace_back(12, NO_LOCREF, 9.0,  3.0, num_SVs_G, &IS_levels[0]);
//  IS_levels.emplace_back(18, NO_LOCREF, 9.0,  6.0, num_SVs_G, &IS_levels[1]);
}