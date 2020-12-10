#include "uncolfalgo.h"

#include <algorithm>

//###################################################################
/**Given a point on the unit sphere this function will first determine
 * the closest Spherical Quadrilateral(SQ) then the closest sub-SQ.
 * After this is determined an index to the point is returned.*/
int UncolFAlgo::GetClosestSSQFromIS(const LDFEIS &IS, const chi_mesh::Vector3 &USP)
{
  //======================================== Compute distances to each SQ
  // Stores the result as a distance-index pair
  typedef std::pair<double,int> DK;
  int num_SQs = IS.base_quadrature.deployed_SQs.size();
  std::vector<DK> d_to_SQs_sqr(num_SQs, DK(1.0e15,-1));

  for (int i=0; i<num_SQs; ++i)
  {
    auto& SQ = IS.base_quadrature.deployed_SQs[i];

    d_to_SQs_sqr[i] = DK((USP - SQ.centroid_xyz).NormSquare(),i);
  }

  //======================================== Sorts the distance-index pairs
  std::sort(d_to_SQs_sqr.begin(),d_to_SQs_sqr.end());

  int ref_SQ_index = d_to_SQs_sqr.front().second;
  auto& ref_SQ = IS.base_quadrature.deployed_SQs[ref_SQ_index];

  //======================================== Compute distances to each sub-SQ
  std::vector<DK> d_to_SSQs_sqr(4,DK(1.0e15,-1));
  for (int i=0; i<4; ++i)
    d_to_SSQs_sqr[i] = DK((USP - ref_SQ.sub_sqr_points[i]).NormSquare(),i);

  //======================================== Sorts
  std::sort(d_to_SSQs_sqr.begin(),d_to_SSQs_sqr.end());

  int ref_SSQ_index = d_to_SSQs_sqr.front().second;

  return ref_SQ_index*4 + ref_SSQ_index;
}