#ifndef GMSH_API_PARALLEL_H
#define GMSH_API_PARALLEL_H

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/tria_description.h>

#ifdef DEAL_II_GMSH_WITH_API
#  include <gmsh.h>
#endif

#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace GMSH
{
#ifdef DEAL_II_GMSH_WITH_API
 // Reads a partitioned GMSH mesh file and creates a fully distributed triangulation in deal.II.
  template <int dim, int spacedim>
  void
  read_partitioned_msh(
    dealii::parallel::fullydistributed::Triangulation<dim, spacedim> &tria,
    const MPI_Comm                                                   &mpi_comm,
    const std::string &file_prefix = "mesh_",
    const std::string &file_suffix = ".msh")
  {
    using namespace dealii;
    const unsigned int nprocs = Utilities::MPI::n_mpi_processes(mpi_comm);
    const unsigned int rank   = Utilities::MPI::this_mpi_process(mpi_comm);

    std::string fname =
      file_prefix + "_" + std::to_string(rank + 1) + file_suffix;

    if (nprocs == 1)
      {
        fname = file_prefix + file_suffix;
        std::ifstream f(fname);
        AssertThrow(f.good(), ExcMessage("Missing mesh file: " + fname));
      }
    else
      {
        // Check that files with names prefix1.suffix to prefixnprocs.suffix
        // exist, and that prefix(nprocs+1).suffix does not exist.
        for (unsigned int i = 1; i <= nprocs; ++i)
          {
            std::string check_fname =
              file_prefix + "_" + std::to_string(i) + file_suffix;
            std::ifstream f(check_fname);
            AssertThrow(f.good(),
                        ExcMessage("Missing mesh file: " + check_fname));
          }

        std::string extra_fname =
          file_prefix + std::to_string(nprocs + 1) + file_suffix;
        std::ifstream f(extra_fname);
        AssertThrow(!f.good(),
                    ExcMessage("Unexpected extra mesh file: " + extra_fname));
      }

    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 0);
    gmsh::open(fname);

    std::map<unsigned long, unsigned int> ghost_map;

    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities);

    for (const auto &e : entities)
    {
        const int entity_dim = e.first;
        const int entity_tag = e.second;

        if (entity_dim == dim)
        {
            std::vector<std::size_t> element_tags;
            std::vector<int> partitions;

            gmsh::model::mesh::getGhostElements(entity_dim, entity_tag, element_tags, partitions);

            for (std::size_t i = 0; i < element_tags.size(); ++i)
            {
                ghost_map[element_tags[i]] = static_cast<unsigned int>(partitions[i] - 1);
            }
        }
    }

    std::vector<std::size_t> node_tags;
    std::vector<double> coords, parametric_coords;
    gmsh::model::mesh::getNodes(node_tags, coords, parametric_coords);

    TriangulationDescription::Description<dim, spacedim> triangulation_description;
    triangulation_description.comm = mpi_comm;

    auto &vertices = triangulation_description.coarse_cell_vertices;
    auto &cells = triangulation_description.coarse_cells;
    auto &coarse_cell_ids = triangulation_description.coarse_cell_index_to_coarse_cell_id;
    auto &cell_infos = triangulation_description.cell_infos;
    cell_infos.resize(1);

    std::map<std::size_t, unsigned int> node_tag_to_index
    vertices.resize(node_tags.size(), Point<spacedim>());
    for (unsigned int i = 0; i < node_tags.size(); ++i)
    {
        node_tag_to_index[node_tags[i]] = i;
        for (unsigned int d = 0; d < spacedim; ++d)
            vertices[i][d] = coords[3 * i + d];
    }

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
                if (element_ids[i].empty())
                    continue;

                const unsigned int n_vertices = element_nodes[i].size() / element_ids[i].size();

                for (unsigned int j = 0; j < element_ids[i].size(); ++j)
                {
                    CellData<dim> cell;
                    cell.material_id = 0;

                    for (unsigned int v = 0; v < n_vertices; ++v)
                    {
                        const std::size_t node_tag = element_nodes[i][j * n_vertices + v];
                        AssertThrow(node_tag_to_index.find(node_tag) != node_tag_to_index.end(),
                                    ExcMessage("Node tag " + std::to_string(node_tag) + " not found in node list!"));
                        cell.vertices[v] = node_tag_to_index[node_tag];
                    }

                    cells.push_back(cell);
                    coarse_cell_ids.push_back(element_ids[i][j]);

                    TriangulationDescription::CellData<dim> cell_info;
                    cell_info.id = CellId(element_ids[i][j], {}).template to_binary<dim>();

                    auto it = ghost_map.find(element_ids[i][j]);
                    if (it != ghost_map.end())
                    {
                        cell_info.subdomain_id = it->second;
                    }
                    else
                    {
                        cell_info.subdomain_id = rank;
                    }

                    cell_info.level_subdomain_id = cell_info.subdomain_id;
                    cell_infos[0].push_back(cell_info);
                }
            }
        }
    }

    triangulation_description.settings = TriangulationDescription::Settings::default_setting;

    tria.create_triangulation(triangulation_description);

    MPI_Barrier(mpi_comm);

    gmsh::clear();
    gmsh::finalize();
  }

#endif
} // namespace GMSH

DEAL_II_NAMESPACE_CLOSE

#endif // GMSH_API_PARALLEL_H
