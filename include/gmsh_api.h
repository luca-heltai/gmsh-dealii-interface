#ifndef GMSH_API_H
#define GMSH_API_H

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/point.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#ifdef DEAL_II_GMSH_WITH_API
#  include <gmsh.h>
#endif

#include <array>
#include <map>
#include <string>
#include <vector>

#include "debug_utils.h" // new


DEAL_II_NAMESPACE_OPEN

namespace GMSH
{
#ifdef DEAL_II_GMSH_WITH_API
  template <int dim, int spacedim>
  void
  read_msh(Triangulation<dim, spacedim> &tria, const std::string &fname)
  {
    // gmsh -> deal.II types
    const std::map<int, std::uint8_t> gmsh_to_dealii_type = {
      {{15, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {7, 5}, {6, 6}, {5, 7}}};

    // Vertex renumbering, by dealii type
    const std::array<std::vector<unsigned int>, 8> gmsh_to_dealii = {
      {{0},
       {{0, 1}},
       {{0, 1, 2}},
       {{0, 1, 3, 2}},
       {{0, 1, 2, 3}},
       {{0, 1, 3, 2, 4}},
       {{0, 1, 2, 3, 4, 5}},
       {{0, 1, 3, 2, 4, 5, 7, 6}}}};

    std::vector<Point<spacedim>>               vertices;
    std::vector<CellData<dim>>                 cells;
    SubCellData                                subcelldata;
    std::map<unsigned int, types::boundary_id> boundary_ids_1d;

    // Track the number of times each vertex is used in 1D. This determines
    // whether or not we can assign a boundary id to a vertex. This is necessary
    // because sometimes gmsh saves internal vertices in the $ELEM list in codim
    // 1 or codim 2.
    std::map<unsigned int, unsigned int> vertex_counts;

    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 0);
    gmsh::open(fname);

    AssertThrow(gmsh::model::getDimension() == dim,
                ExcMessage(
                  "You are trying to read a gmsh file with dimension " +
                  std::to_string(gmsh::model::getDimension()) +
                  " into a grid of dimension " + std::to_string(dim)));

    // Read all nodes, and store them in our vector of vertices. Before we do
    // that, make sure all tags are consecutive
    {
      gmsh::model::mesh::removeDuplicateNodes();
      gmsh::model::mesh::renumberNodes();
      std::vector<std::size_t> node_tags;
      std::vector<double>      coord;
      std::vector<double>      parametricCoord;
      gmsh::model::mesh::getNodes(
        node_tags, coord, parametricCoord, -1, -1, false, false);
      vertices.resize(node_tags.size());
      for (unsigned int i = 0; i < node_tags.size(); ++i)
        {
          // Check that renumbering worked!
          AssertDimension(node_tags[i], i + 1);
          for (unsigned int d = 0; d < spacedim; ++d)
            vertices[i][d] = coord[i * 3 + d];
          if constexpr (running_in_debug_mode())
            {
              // Make sure the embedded dimension is right
              for (unsigned int d = spacedim; d < 3; ++d)
                Assert(std::abs(coord[i * 3 + d]) < 1e-10,
                       ExcMessage(
                         "The grid you are reading contains nodes that are "
                         "nonzero in the coordinate with index " +
                         std::to_string(d) +
                         ", but you are trying to save "
                         "it on a grid embedded in a " +
                         std::to_string(spacedim) + " dimensional space."));
            }
        }
    }

    // Get all the elementary entities in the model, as a vector of (dimension,
    // tag) pairs:
    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities);

    for (const auto &e : entities)
      {
        // Dimension and tag of the entity:
        const int &entity_dim = e.first;
        const int &entity_tag = e.second;

        types::manifold_id manifold_id = numbers::flat_manifold_id;
        types::boundary_id boundary_id = 0;

        // Get the physical tags, to deduce boundary, material, and manifold_id
        std::vector<int> physical_tags;
        gmsh::model::getPhysicalGroupsForEntity(entity_dim,
                                                entity_tag,
                                                physical_tags);

        // Now fill manifold id and boundary or material id
        if (physical_tags.size())
          for (auto physical_tag : physical_tags)
            {
              std::string name;
              gmsh::model::getPhysicalName(entity_dim, physical_tag, name);
              if (!name.empty())
                {
                  // Patterns::Tools::to_value throws an exception, if it can
                  // not convert name to a map from string to int.
                  try
                    {
                      std::map<std::string, int> id_names;
                      Patterns::Tools::to_value(name, id_names);
                      bool found_unrecognized_tag = false;
                      bool found_boundary_id      = false;
                      // If the above did not throw, we keep going, and retrieve
                      // all the information that we know how to translate.
                      for (const auto &it : id_names)
                        {
                          const auto &name = it.first;
                          const auto &id   = it.second;
                          if (entity_dim == dim && name == "MaterialID")
                            {
                              boundary_id = static_cast<types::boundary_id>(id);
                              found_boundary_id = true;
                            }
                          else if (entity_dim < dim && name == "BoundaryID")
                            {
                              boundary_id = static_cast<types::boundary_id>(id);
                              found_boundary_id = true;
                            }
                          else if (name == "ManifoldID")
                            manifold_id = static_cast<types::manifold_id>(id);
                          else
                            // We did not recognize one of the keys. We'll fall
                            // back to setting the boundary id to the physical
                            // tag after reading all strings.
                            found_unrecognized_tag = true;
                        }
                      // If we didn't find a BoundaryID:XX or MaterialID:XX, and
                      // something was found but not recognized, then we set the
                      // material id or boundary using the physical tag
                      // directly.
                      if (found_unrecognized_tag && !found_boundary_id)
                        boundary_id = physical_tag;
                    }
                  catch (...)
                    {
                      // When the above didn't work, we revert to the old
                      // behaviour: the physical tag itself is interpreted
                      // either as a material_id or a boundary_id, and no
                      // manifold id is known
                      boundary_id = physical_tag;
                    }
                }
            }

        // Get the mesh elements for the entity (dim, tag):
        std::vector<int>                      element_types;
        std::vector<std::vector<std::size_t>> element_ids, element_nodes;
        gmsh::model::mesh::getElements(
          element_types, element_ids, element_nodes, entity_dim, entity_tag);

        for (unsigned int i = 0; i < element_types.size(); ++i)
          {
            const auto &type       = gmsh_to_dealii_type.at(element_types[i]);
            const auto  n_vertices = gmsh_to_dealii[type].size();
            const auto &elements   = element_ids[i];
            const auto &nodes      = element_nodes[i];
            for (unsigned int j = 0; j < elements.size(); ++j)
              {
                if (entity_dim == dim)
                  {
                    cells.emplace_back(n_vertices);
                    auto &cell = cells.back();
                    for (unsigned int v = 0; v < n_vertices; ++v)
                      {
                        cell.vertices[v] =
                          nodes[n_vertices * j + gmsh_to_dealii[type][v]] - 1;
                        if (dim == 1)
                          vertex_counts[cell.vertices[v]] += 1u;
                      }
                    cell.manifold_id = manifold_id;
                    cell.material_id = boundary_id;
                  }
                else if (entity_dim == 2)
                  {
                    subcelldata.boundary_quads.emplace_back(n_vertices);
                    auto &face = subcelldata.boundary_quads.back();
                    for (unsigned int v = 0; v < n_vertices; ++v)
                      face.vertices[v] =
                        nodes[n_vertices * j + gmsh_to_dealii[type][v]] - 1;

                    face.manifold_id = manifold_id;
                    face.boundary_id = boundary_id;
                  }
                else if (entity_dim == 1)
                  {
                    subcelldata.boundary_lines.emplace_back(n_vertices);
                    auto &line = subcelldata.boundary_lines.back();
                    for (unsigned int v = 0; v < n_vertices; ++v)
                      line.vertices[v] =
                        nodes[n_vertices * j + gmsh_to_dealii[type][v]] - 1;

                    line.manifold_id = manifold_id;
                    line.boundary_id = boundary_id;
                  }
                else if (entity_dim == 0)
                  {
                    // This should only happen in one dimension.
                    AssertDimension(dim, 1);
                    for (unsigned int j = 0; j < elements.size(); ++j)
                      boundary_ids_1d[nodes[j] - 1] = boundary_id;
                  }
              }
          }
      }

    // apply_grid_fixup_functions() may invalidate the vertex indices (since it
    // will delete duplicated or unused vertices). Get around this by storing
    // Points directly in that case so that the comparisons are valid.
    std::vector<std::pair<Point<spacedim>, types::boundary_id>>
      boundary_id_pairs;
    if (dim == 1)
      for (const auto &pair : vertex_counts)
        if (pair.second == 1u)
          boundary_id_pairs.emplace_back(vertices[pair.first],
                                         boundary_ids_1d[pair.first]);

    // apply_grid_fixup_functions(vertices, cells, subcelldata);
    tria.create_triangulation(vertices, cells, subcelldata);

    // // in 1d, we also have to attach boundary ids to vertices, which does not
    // // currently work through the call above.
    // if (dim == 1)
    //   assign_1d_boundary_ids(boundary_id_pairs, tria);

    gmsh::clear();
    gmsh::finalize();
  }
#endif
} // namespace GMSH

DEAL_II_NAMESPACE_CLOSE

//

#endif // GMSH_API_H
