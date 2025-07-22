#ifndef GMSH_API_SERIALIZATION_H
#define GMSH_API_SERIALIZATION_H

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
#include <map>
#include <sstream>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace GMSH
{
#ifdef DEAL_II_GMSH_WITH_API

  // Standalone function, not nested!
  inline std::map<unsigned long, unsigned int>
  parse_ghost_elements(const std::string &filename)
  {
    std::map<unsigned long, unsigned int> ghost_map;
    std::ifstream                         in(filename);
    std::string                           line;
    bool                                  in_ghost = false;
    while (std::getline(in, line))
      {
        if (line.find("$GhostElements") != std::string::npos)
          {
            in_ghost = true;
            continue;
          }
        if (in_ghost && line.find("$EndGhostElements") != std::string::npos)
          break;
        if (in_ghost && !line.empty() && std::isdigit(line[0]))
          {
            std::istringstream iss(line);
            unsigned long      cell_id;
            int                dim, owner, ghost_on;
            if (iss >> cell_id >> dim >> owner >> ghost_on)
              ghost_map[cell_id] = owner - 1; // GMSH 1-based to deal.II 0-based
          }
      }
    return ghost_map;
  }

  template <int dim, int spacedim>
  void
  read_parallel_msh(
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

    // Parse ghost elements before GMSH API
    std::map<unsigned long, unsigned int> ghost_map =
      parse_ghost_elements(fname);

    // Variables to hold mesh data
    std::vector<std::size_t>             node_tags;
    std::vector<double>                  coords, parametric_coords;
    std::vector<std::pair<int, int>>     entities;
    std::map<std::size_t, unsigned int>  node_tag_to_index;
    std::vector<dealii::Point<spacedim>> vertices;
    std::vector<dealii::CellData<dim>>   cells;
    std::vector<unsigned long>           coarse_cell_ids;
    std::vector<std::vector<dealii::TriangulationDescription::CellData<dim>>>
      cell_infos(1);

    // --- SERIALIZE GMSH API USAGE ---
    for (unsigned int r = 0; r < nprocs; ++r)
      {
        if (rank == r)
          {
            std::cout << "Rank " << rank << " initializing GMSH" << std::endl;
            gmsh::initialize();
            gmsh::option::setNumber("General.Verbosity", 0);
            std::cout << "Rank " << rank << " opening mesh: " << fname
                      << std::endl;
            gmsh::open(fname);

            AssertThrow(gmsh::model::getDimension() == dim,
                        ExcMessage("Dimension mismatch: Gmsh file dimension " +
                                   std::to_string(gmsh::model::getDimension()) +
                                   " vs triangulation dimension " +
                                   std::to_string(dim)));

            // Retrieve nodes
            gmsh::model::mesh::getNodes(node_tags, coords, parametric_coords);

            // Map node tags to indices and fill vertices
            node_tag_to_index.clear();
            vertices.resize(node_tags.size(), dealii::Point<spacedim>());
            for (unsigned int i = 0; i < node_tags.size(); ++i)
              {
                node_tag_to_index[node_tags[i]] = i;
                for (unsigned int d = 0; d < spacedim; ++d)
                  vertices[i][d] = coords[3 * i + d];
              }

            // Retrieve elements
            gmsh::model::getEntities(entities);

            cells.clear();
            coarse_cell_ids.clear();
            cell_infos[0].clear();

            for (const auto &e : entities)
              {
                const int entity_dim = e.first;
                const int entity_tag = e.second;

                if (entity_dim == dim)
                  {
                    std::vector<int>                      element_types;
                    std::vector<std::vector<std::size_t>> element_ids,
                      element_nodes;
                    gmsh::model::mesh::getElements(element_types,
                                                   element_ids,
                                                   element_nodes,
                                                   entity_dim,
                                                   entity_tag);

                    for (unsigned int i = 0; i < element_types.size(); ++i)
                      {
                        if (element_ids[i].empty())
                          continue;

                        const unsigned int n_vertices =
                          element_nodes[i].size() / element_ids[i].size();

                        for (unsigned int j = 0; j < element_ids[i].size(); ++j)
                          {
                            dealii::CellData<dim> cell;
                            cell.material_id = 0;

                            for (unsigned int v = 0; v < n_vertices; ++v)
                              {
                                const std::size_t node_tag =
                                  element_nodes[i][j * n_vertices + v];
                                AssertThrow(node_tag_to_index.find(node_tag) !=
                                              node_tag_to_index.end(),
                                            dealii::ExcMessage(
                                              "Node tag " +
                                              std::to_string(node_tag) +
                                              " not found in node list!"));
                                cell.vertices[v] = node_tag_to_index[node_tag];
                              }

                            cells.push_back(cell);
                            coarse_cell_ids.push_back(element_ids[i][j]);

                            // Assign subdomain_id: if ghost, use owner; else,
                            // local rank
                            dealii::TriangulationDescription::CellData<dim>
                              cell_info;
                            cell_info.id = dealii::CellId(element_ids[i][j], {})
                                             .template to_binary<dim>();

                            auto it = ghost_map.find(element_ids[i][j]);
                            if (it != ghost_map.end())
                              cell_info.subdomain_id =
                                it->second; // Owner from ghost map
                            else
                              cell_info.subdomain_id = rank; // Locally owned

                            cell_info.level_subdomain_id =
                              cell_info.subdomain_id;
                            cell_infos[0].push_back(cell_info);
                          }
                      }
                  }
              }

            gmsh::clear();
            gmsh::finalize();
            std::cout << "Rank " << rank << " finalized GMSH" << std::endl;
          }
        MPI_Barrier(mpi_comm);
      }

    // Now, all ranks have filled vertices, cells, coarse_cell_ids, cell_infos

    // Prepare triangulation description
    dealii::TriangulationDescription::Description<dim, spacedim>
      triangulation_description;
    triangulation_description.comm                 = mpi_comm;
    triangulation_description.coarse_cell_vertices = vertices;
    triangulation_description.coarse_cells         = cells;
    triangulation_description.coarse_cell_index_to_coarse_cell_id =
      coarse_cell_ids;
    triangulation_description.cell_infos = cell_infos;
    triangulation_description.settings =
      dealii::TriangulationDescription::Settings::default_setting;

    tria.create_triangulation(triangulation_description);
  }
#endif
} // namespace GMSH

DEAL_II_NAMESPACE_CLOSE

#endif // GMSH_API_SERIALIZATION_H
