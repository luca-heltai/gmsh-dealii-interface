#ifndef GMSH_API_PARALLEL_H
#define GMSH_API_PARALLEL_H

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/point.h>
#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/grid/tria_description.h>
#include <iostream>
#include <csignal>
#include <cmath>

#ifdef DEAL_II_GMSH_WITH_API
#  include <gmsh.h>
#endif

#include <array>
#include <map>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace GMSH
{
#ifdef DEAL_II_GMSH_WITH_API
  /**
   * @brief Read partitioned Gmsh mesh files (one per MPI rank)
   *        and populate a parallel::fullydistributed::Triangulation.
   */
  template <int dim, int spacedim>
  void read_parallel_msh(
    dealii::parallel::fullydistributed::Triangulation<dim, spacedim> &tria,
    const MPI_Comm                                                   &mpi_comm,
    const std::string &file_prefix = "mesh_",
    const std::string &file_suffix = ".msh")
  {
    using namespace dealii;

    const unsigned int rank = Utilities::MPI::this_mpi_process(mpi_comm);
    const std::string fname = file_prefix + std::to_string(rank) + file_suffix;

    // Initialize Gmsh API and suppress verbosity
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 1); // Set verbosity level for debugging
    gmsh::open(fname);

    // Check that mesh dimension matches expected dimension
    AssertThrow(gmsh::model::getDimension() == dim,
                ExcMessage("Gmsh file dimension does not match expected dim."));

    // Get all node tags and coordinates from the mesh
    std::vector<std::size_t> node_tags;
    std::vector<double> coords, parametric_coords;
    gmsh::model::mesh::getNodes(node_tags, coords, parametric_coords);

    // Debug: Check node_tags and coords size
    AssertThrow(node_tags.size() * 3 == coords.size(),
                ExcMessage("Mismatch between node tags and coordinates size"));

    // Debug: Check for non-finite coordinates
    for (unsigned int i = 0; i < node_tags.size(); ++i)
    {
      for (unsigned int d = 0; d < spacedim; ++d)
      {
        double coord = coords[i * 3 + d];
        if (!std::isfinite(coord))
        {
          std::cerr << "Rank " << rank << ": Non-finite coordinate at node "
                    << node_tags[i] << " dim " << d << ": " << coord << std::endl;
          AssertThrow(false, ExcMessage("Non-finite coordinate in mesh"));
        }
      }
    }

    // Prepare the triangulation description
    TriangulationDescription::Description<dim, spacedim> triangulation_description;
    triangulation_description.comm = mpi_comm;

    auto &vertices = triangulation_description.coarse_cell_vertices;
    auto &cells = triangulation_description.coarse_cells;
    auto &coarse_cell_ids = triangulation_description.coarse_cell_index_to_coarse_cell_id;

    // Map node tag → local index, and fill vertex coordinates
    std::map<std::size_t, unsigned int> node_tag_to_index;
    vertices.resize(node_tags.size());

    for (unsigned int i = 0; i < node_tags.size(); ++i)
    {
      node_tag_to_index[node_tags[i]] = i;
      for (unsigned int d = 0; d < spacedim; ++d)
        vertices[i][d] = coords[i * 3 + d];  // Gmsh stores 3D coords
    }

    // Mapping of Gmsh element type → deal.II cell type (numbered 0..7)
    const std::map<int, std::uint8_t> gmsh_to_dealii_type = {
      {15, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {7, 5}, {6, 6}, {5, 7}
    };

    // Vertex reordering: Gmsh vertex numbering → deal.II vertex ordering
    const std::array<std::vector<unsigned int>, 8> gmsh_to_dealii = {{
      {0},                     // 0-node
      {0, 1},                  // line
      {0, 1, 2},               // triangle
      {0, 1, 3, 2},            // quadrilateral
      {0, 1, 2, 3},            // tetrahedron
      {0, 1, 3, 2, 4},         // pyramid
      {0, 1, 2, 3, 4, 5},      // prism
      {0, 1, 3, 2, 4, 5, 7, 6} // hexahedron
    }};

    // Loop over all topological entities in the mesh
    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities);

    for (const auto &[entity_dim, entity_tag] : entities)
    {
      // Skip anything not of top dimension (volume cells only)
      if (entity_dim != dim)
        continue;

      std::vector<int> element_types;
      std::vector<std::vector<std::size_t>> element_ids, element_nodes;
      gmsh::model::mesh::getElements(element_types, element_ids, element_nodes,
                                     entity_dim, entity_tag);

      for (unsigned int t = 0; t < element_types.size(); ++t)
      {
        if (element_ids[t].empty() || element_nodes[t].empty())
        {
            std::cerr << "Rank: " << rank << " - Warning: Empty element_ids or element_nodes for entity " << entity_tag << std::endl;
            continue;
        }

        const int gmsh_type = element_types[t];

        auto it = gmsh_to_dealii_type.find(gmsh_type);
        AssertThrow(it != gmsh_to_dealii_type.end(),
                    ExcMessage("Unsupported Gmsh element type: " + std::to_string(gmsh_type)));
        const unsigned int deal_type = it->second;
        const unsigned int n_vertices = gmsh_to_dealii[deal_type].size();

        AssertThrow(n_vertices > 0, ExcMessage("Invalid number of vertices for element type."));

        const auto &ids = element_ids[t];
        const auto &nodes = element_nodes[t];

        for (unsigned int j = 0; j < ids.size(); ++j)
        {
          typename TriangulationDescription::Description<dim, spacedim>::CoarseCell cell;
          cell.vertices.resize(n_vertices);

          for (unsigned int v = 0; v < n_vertices; ++v)
          {
            const std::size_t gmsh_node_tag = nodes[j * n_vertices + gmsh_to_dealii[deal_type][v]];
            if (node_tag_to_index.find(gmsh_node_tag) == node_tag_to_index.end())
            {
              std::cerr << "Rank " << rank << ": Invalid node tag in element: " << gmsh_node_tag << std::endl;
              AssertThrow(false, ExcMessage("Invalid node tag in Gmsh mesh."));
            }
            cell.vertices[v] = node_tag_to_index.at(gmsh_node_tag);
          }

          // Assign dummy CellId (it will be fixed later internally)
          cell.cell_id = CellId::invalid_cell_id();

          cells.push_back(std::move(cell));
          coarse_cell_ids.push_back(ids[j]); // unique coarse cell identifier
        }
      }
    }

    // Debug: Check cell vertex indices are valid
    for (const auto &cell : cells)
    {
      for (auto v : cell.vertices)
      {
        if (v >= vertices.size())
        {
          std::cerr << "Rank " << rank << ": Cell vertex index out of range: " << v 
                    << " (vertices.size() = " << vertices.size() << ")" << std::endl;
          AssertThrow(false, ExcMessage("Cell vertex index out of range"));
        }
      }
    }

    // Now generate the actual parallel triangulation from the description
    AssertThrow(!cells.empty(), ExcMessage("No cells found in triangulation description."));
    AssertThrow(!vertices.empty(), ExcMessage("No vertices found in triangulation description."));

    // Add this before create_triangulation()
std::cout << "Rank " << rank << ": "
<< "Vertices: " << vertices.size() << ", "
<< "Cells: " << cells.size() << ", "
<< "Coarse cell IDs: " << coarse_cell_ids.size() << std::endl;

// Verify description consistency
 AssertThrow(cells.size() == coarse_cell_ids.size(),
  ExcMessage("Cell count doesn't match ID count"));

    tria.create_triangulation(triangulation_description);

    gmsh::clear();
    gmsh::finalize();
  }
#endif
} // namespace GMSH

DEAL_II_NAMESPACE_CLOSE

#endif // GMSH_API_PARALLEL_H
