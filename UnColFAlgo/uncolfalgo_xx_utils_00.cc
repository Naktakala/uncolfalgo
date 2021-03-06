#include "uncolfalgo.h"
#include "ChiMesh/Cell/cell_polygon.h"
#include "ChiMesh/Cell/cell_polyhedron.h"

//###################################################################
/**Determines the centroid of a list of vertices.*/
chi_mesh::Vector3 UncolFAlgo::GetCentroidFromList(const VertList &vert_list)
{
  chi_mesh::Vector3 centroid;

  for (const auto& vert : vert_list)
    centroid += vert;
  centroid /= vert_list.size();

  return centroid;
}

//###################################################################
/**Recursively subdivides a triangles.*/
UncolFAlgo::VertCollection UncolFAlgo::
  SubdivideTriangle(const VertCollection &current_triangle_coll,
                    const int desired_level)
{
  int current_level = log(current_triangle_coll.size())/log(4.0);
  if (current_level >= desired_level)
    return current_triangle_coll;

  VertCollection new_collection;
  for (const auto& triangle : current_triangle_coll)
  {
    auto& v0 = triangle[0];
    auto& v1 = triangle[1];
    auto& v2 = triangle[2];

    auto v01 = 0.5*(v0+v1);
    auto v12 = 0.5*(v1+v2);
    auto v20 = 0.5*(v2+v0);

    new_collection.push_back({v0 ,v01,v20});
    new_collection.push_back({v01,v1 ,v12});
    new_collection.push_back({v12,v2 ,v20});
    new_collection.push_back({v20,v01,v12});
  }

  int new_level = log(new_collection.size())/log(4.0);

  if (new_level >= desired_level)
    return new_collection;
  else
    return SubdivideTriangle(new_collection,desired_level);
}

//###################################################################
/**Spawns a triangle from a set of vertex ids.*/
chi_mesh::Cell* UncolFAlgo::SpawnTriangle(const VertList& verts,
                                          chi_mesh::MeshContinuum* grid)
{
  auto khat = chi_mesh::Vector3(0.0,0.0,1.0);

  const int prev_last_vid = grid->vertices.size()-1;

  grid->vertices.push_back(new chi_mesh::Vertex(verts[0]));
  grid->vertices.push_back(new chi_mesh::Vertex(verts[1]));
  grid->vertices.push_back(new chi_mesh::Vertex(verts[2]));

  std::vector<int> vids = {prev_last_vid+1,
                           prev_last_vid+2,
                           prev_last_vid+3};

  auto& v0 = verts[0];
  auto& v1 = verts[1];
  auto& v2 = verts[2];

  auto new_cell = new chi_mesh::CellPolygon();

  new_cell->vertex_ids = {vids[0],vids[1],vids[2]};
  new_cell->centroid = GetCentroidFromList({v0,v1,v2});

  new_cell->faces.resize(3);
  new_cell->faces[0].vertex_ids = {vids[0],vids[1]};
  new_cell->faces[1].vertex_ids = {vids[1],vids[2]};
  new_cell->faces[2].vertex_ids = {vids[2],vids[0]};

  new_cell->faces[0].centroid = 0.5*(v0+v1);
  new_cell->faces[1].centroid = 0.5*(v1+v2);
  new_cell->faces[2].centroid = 0.5*(v2+v0);

  new_cell->faces[0].normal = (v1-v0).Cross(khat).Normalized();
  new_cell->faces[1].normal = (v2-v1).Cross(khat).Normalized();
  new_cell->faces[2].normal = (v0-v2).Cross(khat).Normalized();

  return new_cell;
}

//###################################################################
/**Recursively subdivides a triangles.*/
UncolFAlgo::VertCollection UncolFAlgo::
SubDivideTetrahedron(const VertCollection &current_tet_coll,
                     const int desired_level)
{
  int current_level = log(current_tet_coll.size())/log(12.0);
  if (current_level >= desired_level)
    return current_tet_coll;

  VertCollection new_collection;
  for (const auto& tet : current_tet_coll)
  {
    VertCollection faces = {{tet[0],tet[1],tet[3]},
                            {tet[0],tet[3],tet[2]},
                            {tet[0],tet[2],tet[1]},
                            {tet[1],tet[2],tet[3]}};

    auto cc = GetCentroidFromList(tet);

    for (const auto& face : faces)
    {
      auto fc = GetCentroidFromList(face);

      auto& v0 = face[0];
      auto& v1 = face[1];
      auto& v2 = face[2];

      new_collection.push_back({cc,fc,v0,v1});
      new_collection.push_back({cc,fc,v1,v2});
      new_collection.push_back({cc,fc,v2,v0});
    }

  }

  int new_level = log(new_collection.size())/log(12.0);

  if (new_level >= desired_level)
    return new_collection;
  else
    return SubdivideTriangle(new_collection,desired_level);
}

//###################################################################
/**Spawns a triangle from a set of vertex ids.*/
chi_mesh::Cell* UncolFAlgo::SpawnTetrahedron(const VertList& verts,
                                             chi_mesh::MeshContinuum* grid)
{
  const int prev_last_vid = grid->vertices.size()-1;

  grid->vertices.push_back(new chi_mesh::Vertex(verts[0]));
  grid->vertices.push_back(new chi_mesh::Vertex(verts[1]));
  grid->vertices.push_back(new chi_mesh::Vertex(verts[2]));
  grid->vertices.push_back(new chi_mesh::Vertex(verts[3]));

  std::vector<int> vids = {prev_last_vid+1,
                           prev_last_vid+2,
                           prev_last_vid+3,
                           prev_last_vid+4};

  auto& v0 = verts[0];
  auto& v1 = verts[1];
  auto& v2 = verts[2];
  auto& v3 = verts[3];

  auto new_cell = new chi_mesh::CellPolyhedron();

  new_cell->vertex_ids = {vids[0],vids[1],vids[2],vids[3]};
  new_cell->centroid = GetCentroidFromList({v0,v1,v2,v3});

  new_cell->faces.resize(4);
  new_cell->faces[0].vertex_ids = {vids[0],vids[2],vids[1]};
  new_cell->faces[1].vertex_ids = {vids[2],vids[3],vids[1]};
  new_cell->faces[2].vertex_ids = {vids[3],vids[0],vids[1]};
  new_cell->faces[3].vertex_ids = {vids[0],vids[3],vids[2]};

  new_cell->faces[0].centroid = GetCentroidFromList({v0,v2,v1});
  new_cell->faces[1].centroid = GetCentroidFromList({v2,v3,v1});
  new_cell->faces[2].centroid = GetCentroidFromList({v3,v0,v1});
  new_cell->faces[3].centroid = GetCentroidFromList({v0,v3,v2});

  new_cell->faces[0].normal = (v2-v0).Cross(v1-v0).Normalized();
  new_cell->faces[1].normal = (v3-v2).Cross(v1-v2).Normalized();
  new_cell->faces[2].normal = (v0-v3).Cross(v1-v3).Normalized();
  new_cell->faces[3].normal = (v3-v0).Cross(v2-v0).Normalized();

  return new_cell;
}