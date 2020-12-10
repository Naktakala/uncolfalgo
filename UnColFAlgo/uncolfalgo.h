#ifndef _uncolfalgo_h
#define _uncolfalgo_h

#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/Cell/cell.h"

#include "ChiMath/SpatialDiscretization/PiecewiseLinear/pwl.h"
#include "ChiMath/UnknownManager/unknown_manager.h"
#include "ChiMath/Quadratures/quadrature_gausslegendre.h"

#include "ChiMath/Quadratures/SLDFESQ/sldfe_sq.h"

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

  //xx_02
  struct LDFEIS;
  int GetClosestSSQFromIS(const LDFEIS& IS,const chi_mesh::Vector3& USP);


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
/**LDFE quadrature refinitement structure.*/
struct LDFELocalRefinement
{
  const chi_mesh::Vector3 ref_dir;
  const double            cone_size;
  const bool              dir_as_plane_normal;

  LDFELocalRefinement(const chi_mesh::Vector3& in_ref_dir,
                      double in_cone_size,
                      bool in_dir_as_plane_normal ) :
                      ref_dir(in_ref_dir),
                      cone_size(in_cone_size),
                      dir_as_plane_normal(in_dir_as_plane_normal)
  { }
};
//###################################################################
/**LDFE Interpolation Surface data structure.*/
struct LDFEIS
{
  friend int UncolFAlgo::GetClosestSSQFromIS(const LDFEIS &IS,
                                             const chi_mesh::Vector3 &USP);
private:
  chi_math::SimplifiedLDFESQ::Quadrature base_quadrature;
  const int                              initial_refinement;
  const std::vector<LDFELocalRefinement> local_refinements;
public:
  const double                           radius_outer_sqr;
  const double                           radius_inner_sqr;
  const LDFEIS*                          parent;

public:
  std::vector<chi_mesh::Vector3>                 cone_directions;
  std::vector<double>                            cone_solid_angles;
  std::vector<std::vector<double>>               SPig_taus;
  std::vector<int>                               parent_SSQ_mapping;
  LDFEIS(int in_initial_refinement,
         std::vector<LDFELocalRefinement> in_local_refinements,
         double in_radius_outer,
         double in_radius_inner,
         int num_SVs_times_g,
         LDFEIS* in_parent) :
    initial_refinement(in_initial_refinement),
    local_refinements(std::move(in_local_refinements)),
    radius_outer_sqr(in_radius_outer*in_radius_outer),
    radius_inner_sqr(
      (in_radius_inner<0.0) ? -1.0e-8 : in_radius_inner*in_radius_inner),
    parent(in_parent)
  {
    base_quadrature.GenerateInitialRefinement(initial_refinement);
    for (auto& loc_ref : local_refinements)
      base_quadrature.LocallyRefine(loc_ref.ref_dir,
                                    loc_ref.cone_size,
                                    loc_ref.dir_as_plane_normal);

    cone_directions     = base_quadrature.omegas;
    cone_solid_angles   = base_quadrature.weights;

    SPig_taus.resize(num_SVs_times_g,VecDbl(cone_directions.size(), 0.0));
    if (in_parent != nullptr)
    {
      int num_dirs = cone_directions.size();
      parent_SSQ_mapping.resize(num_dirs,-1);

      for (int i=0; i<num_dirs; ++i)
        parent_SSQ_mapping[i] =
          UncolFAlgo::GetClosestSSQFromIS(*parent, cone_directions[i]);
    }
  }
};

//###################################################################
/**Ray data structure.*/
struct Ray
{
  chi_mesh::Vector3 omega;
  const double      solid_angle=0.0;
  int               current_cell=-1;
  chi_mesh::Vector3 current_position;
  double            tau=0.0;
  const int         IS_level;
  const int         ISP_index;
  const int         SV_index;

  chi_mesh::Vector3 exit_point;

  Ray(const chi_mesh::Vector3& direction,
      double in_solid_angle,
      int current_cell_global_id,
      const chi_mesh::Vector3& in_current_position,
      int in_IS_level,
      int in_ISP_index,
      int in_SV_index) :
    omega(direction),
    solid_angle(in_solid_angle),
    current_cell(current_cell_global_id),
    current_position(in_current_position),
    IS_level(in_IS_level),
    ISP_index(in_ISP_index),
    SV_index(in_SV_index)
  {}
};

//###################################################################
/**Solver class.*/
class Solver
{
private:
  const int                        num_groups;
  const int                        num_moments;
  chi_mesh::MeshContinuum*         grid;
  GeometryType                     geometry_type = GeometryType::NO_GEOMETRY_SET;
  std::vector<chi_mesh::CellFace>  IS_faces;
  VertCollection                   ISs;
  chi_math::SimplifiedLDFESQ::Quadrature LDFEISs_temp;
  std::vector<LDFEIS>              IS_levels;
  std::vector<SourceVolume>        SVs;

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
  void BuildISs();
  void BuildSVs();
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