#include "uncolfalgo.h"

//###################################################################
/**Computes the ray solid angle using spherical excess.*/
double UncolFAlgo::ComputeSolidAngle(
  const chi_mesh::Vector3& SP,
  const VertList& surface_vertices,
  GeometryType geometry_type)
{
  double solid_angle = 0.0;
  if (geometry_type == GeometryType::ONED_SLAB)
    solid_angle = ComputeRaySolidAngle1D(SP,surface_vertices);
  else if (geometry_type == GeometryType::TWOD_CARTESIAN)
    solid_angle = ComputeRaySolidAngle2D(SP,surface_vertices);
  else if (geometry_type == GeometryType::THREED_CARTESIAN)
    solid_angle =  ComputeRaySolidAngle3D(SP,surface_vertices);

  return solid_angle;
}

//###################################################################
/**Computes the ray solid angle using spherical excess.*/
double UncolFAlgo::ComputeRaySolidAngle1D(
  const chi_mesh::Vector3& SP,
  const VertList& surface_vertices)
{
  return 2.0*M_PI;
}

//###################################################################
/**Computes the ray solid angle using spherical excess.*/
double UncolFAlgo::ComputeRaySolidAngle2D(
  const chi_mesh::Vector3& SP,
  const VertList& surface_vertices)
{
  auto& v0 = surface_vertices[0];
  auto& v1 = surface_vertices[1];

  auto vc0 = (v0 - SP).Normalized();
  auto vc1 = (v1 - SP).Normalized();

  double mu = std::max(-1.0,std::min(1.0,vc0.Dot(vc1)));

  if (vc0.Cross(vc1).z > 0.0)
    return  std::fabs(acos(mu)*2.0);
  else
    return -std::fabs(acos(mu)*2.0);
}

//###################################################################
/**Computes the ray solid angle using spherical excess.*/
double UncolFAlgo::ComputeRaySolidAngle3D(
  const chi_mesh::Vector3& SP,
  const VertList& surface_vertices)
{
  //====================================== Lambda for spherical excess
  auto GetSphericalExcess = [](const chi_mesh::Vector3& vA,
                               const chi_mesh::Vector3& vB,
                               const chi_mesh::Vector3& vC)
  {
    const auto& n = vA;

    auto vAB = vB - vA;
    auto vAC = vC - vA;

    auto tAB = vAB.Cross(n);
    auto tAC = vAC.Cross(n);

    auto bAB = n.Cross(tAB).Normalized();
    auto bAC = n.Cross(tAC).Normalized();

    double mu = std::max(-1.0,std::fmin(1.0,bAB.Dot(bAC)));

    return acos(mu);
  };

  double excess = 0.0;
  excess += GetSphericalExcess(surface_vertices[0]-SP,
                               surface_vertices[1]-SP,
                               surface_vertices[2]-SP);
  excess += GetSphericalExcess(surface_vertices[1]-SP,
                               surface_vertices[2]-SP,
                               surface_vertices[0]-SP);
  excess += GetSphericalExcess(surface_vertices[2]-SP,
                               surface_vertices[0]-SP,
                               surface_vertices[1]-SP);

  return excess - M_PI;
}

