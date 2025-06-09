#ifndef GMSH_API_PARALLEL_H
#define GMSH_API_PARALLEL_H

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/point.h>
#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/tria.h>

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
  template <int dim, int spacedim>
  void read_parallel_msh(
    dealii::parallel::fullydistributed::Triangulation<dim, spacedim> &tria,
    const MPI_Comm &  mpi_comm,
    const std::string &  file_prefix = "mesh_",
    const std::string & file_suffix = ".msh")
  {
    using namespace dealii;

    const unsigned int rank = Utilities::MPI::this_mpi_process(mpi_comm);
    const std::string  fname = file_prefix + std::to_string(rank) + file_suffix;

// Initialize GMSH and read the mesh file
    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 0);
    gmsh::open(fname);
    AssertThrow(gmsh::model::getDimension() == dim,
                ExcMessage(
                  "You are trying to read a gmsh file with dimension " +
                  std::to_string(gmsh::model::getDimension()) +
                  " into a grid of dimension " + std::to_string(dim)));
    
    // Read nodes
    std::vector<std::size_t> node_tags;
    std::vector<double> coord;
    gmsh::model::mesh::getNodes(node_tags, coord);

    std::vector<Point<spacedim>> vertices(node_tags.size());
    for (unsigned int i = 0; i < node_tags.size(); ++i)
      for (unsigned int d = 0; d < spacedim; ++d)
        vertices[i][d] = coord[i * 3 + d];

        // Read Elements
        std::vector<std::pair<int, int>> entities;
        gmsh::model::getEntities(entities);

        std::vector<CellData<dim>> cells;
    SubCellData subcelldata;

    for (const auto &e : entities)
    {
      const int entity_dim = e.first;
      const int entity_tag = e.second;

      if (entity_dim == dim)
      {
        std::vector<int> element_types;
        std::vector<std::vector<std::size_t>> element_ids, element_nodes;
        gmsh::model::mesh::getElements(element_types, element_ids, element_nodes, entity_dim, entity_tag);

        for (unsigned int i = 0; i < element_types.size(); ++i)
        {
          const unsigned int n_vertices = element_nodes[i].size() / element_ids[i].size();
          for (unsigned int j = 0; j < element_ids[i].size(); ++j)
          {
            CellData<dim> cell;
            cell.material_id = 0; // Set material ID if needed
            for (unsigned int v = 0; v < n_vertices; ++v)
              cell.vertices[v] = element_nodes[i][j * n_vertices + v] - 1;
            cells.push_back(cell);
          }
        }
      }
    }

    // Partition the mesh
    std::vector<types::subdomain_id> partitioning(cells.size(), rank);

    // Create the distributed triangulation
    tria.create_triangulation(vertices, cells, subcelldata, partitioning);

    // Finalize Gmsh
    gmsh::clear();
    gmsh::finalize();
  }
}

