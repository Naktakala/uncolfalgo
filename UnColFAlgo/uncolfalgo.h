#ifndef _uncolfalgo_h
#define _uncolfalgo_h

#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/Cell/cell.h"

#include "ChiMath/SpatialDiscretization/PiecewiseLinear/pwl.h"
#include "ChiMath/UnknownManager/unknown_manager.h"
#include "ChiMath/Quadratures/quadrature_gausslegendre.h"

namespace UncolFAlgo
{
  typedef std::vector<chi_mesh::Vertex> VertList;
  typedef std::vector<VertList> VertCollection;

  enum class GeometryType
  {
    NO_GEOMETRY_SET  = 0,
    ONED_SLAB        = 1,
    ONED_SPHERICAL   = 2,
    TWOD_CARTESIAN   = 3,
    TWOD_CYLINDRICAL = 4,
    THREED_CARTESIAN = 5
  };

  //xx_00
  chi_mesh::Vector3 GetCentroidFromList(const VertList& vert_list);

  VertCollection SubdivideTriangle(const VertCollection& current_triangle_coll,
                                   const int desired_level);
  chi_mesh::Cell* SpawnTriangle(const VertList& verts,
                                chi_mesh::MeshContinuum* grid);

  VertCollection SubDivideTetrahedron(const VertCollection& current_tet_coll,
                                      const int desired_level);
  chi_mesh::Cell* SpawnTetrahedron(const VertList& verts,
                                   chi_mesh::MeshContinuum* grid);



  //xx_01
  double ComputeSolidAngle(
    const chi_mesh::Vector3& SP,
    const VertList& surface_vertices,
    GeometryType geometry_type);

  double ComputeRaySolidAngle1D(
    const chi_mesh::Vector3& SP,
    const VertList& surface_vertices);

  double ComputeRaySolidAngle2D(
    const chi_mesh::Vector3& SP,
    const VertList& surface_vertices);

  double ComputeRaySolidAngle3D(
    const chi_mesh::Vector3& SP,
    const VertList& surface_vertices);


//###################################################################
/**Source point data structure.*/
struct SourceVolume
{
  const chi_mesh::Vector3 SP; ///< Source Point
  const double magnitude;
  const int owning_cell_global_id;

  const VertCollection faces;

  const double sig_t;

  chi_mesh::Cell* ref_cell = nullptr;

  SourceVolume(const chi_mesh::Vector3& in_source_point,
               const double in_magnitude,
               const int in_owning_cell_global_id,
               const VertCollection& in_faces,
               double in_sig_t) :
    SP(in_source_point),
    magnitude(in_magnitude),
    owning_cell_global_id(in_owning_cell_global_id),
    faces(in_faces),
    sig_t(in_sig_t)
  {
  }

  SourceVolume(const SourceVolume& other) :
    SP(other.SP),
    magnitude(other.magnitude),
    owning_cell_global_id(other.owning_cell_global_id),
    faces(other.faces),
    sig_t(other.sig_t),
    ref_cell(other.ref_cell)
  {
  }
};


//###################################################################
/**Ray data structure.*/
struct Ray
{
  chi_mesh::Vector3 omega;
  double            solid_angle=0.0;
  int               current_cell=-1;
  chi_mesh::Vector3 current_position;
  double            tau=0.0;
  const SourceVolume&     SV;

  chi_mesh::Vector3 exit_point;

  Ray(const SourceVolume& in_SV,
      const chi_mesh::Vector3& direction,
      double in_solid_angle,
      int current_cell_global_id,
      const chi_mesh::Vector3& in_current_position) :
    SV(in_SV),
    omega(direction),
    solid_angle(in_solid_angle),
    current_cell(current_cell_global_id),
    current_position(in_current_position)
  {

  }
};

//###################################################################
/**Solver class.*/
class Solver
{
private:
  const int num_groups;
  const int num_moments;
  chi_mesh::MeshContinuum* grid;
  GeometryType geometry_type = GeometryType::NO_GEOMETRY_SET;
  std::vector<chi_mesh::CellFace> IS_faces;
  VertCollection              ISs;
  std::vector<SourceVolume>   SVs;

  SpatialDiscretization_PWL* pwl;
  chi_math::UnknownManager uk_man;

  unsigned long local_dof_count=0;
  unsigned long globl_dof_count=0;

  chi_math::QuadratureGaussLegendre legendre_quadrature;
  std::vector<double> segment_lengths;
  std::vector<double> shape_values_i;
  std::vector<double> shape_values_f;

  std::vector<double> phi_old_local;
  std::vector<double> phi_old_local_src;
  std::vector<MatDbl> Ak;
  std::vector<double> phi_new_local;

public:
  Solver(const int in_num_groups, const int in_num_moments) :
  num_groups(in_num_groups),
  num_moments(in_num_moments)
  {
    legendre_quadrature.Initialize(3);
  }
  //01
  void Initialize();
  //02
  void Execute();
  //02a
  bool RayTrace(Ray& ray);
  //02b
  void ProcessTrackOutsideSV(Ray& ray,
                             const chi_mesh::Vector3& ray_destination,
                             chi_mesh::Cell& cell,
                             double sig_t);
  //02c
  void ProcessTrackInsideSV(Ray& ray,
                            const chi_mesh::Vector3& ray_destination,
                            chi_mesh::Cell& cell,
                            double sig_t);

  //02d
  void ComputePWLDTransformations();
};

}//namespace

#endif