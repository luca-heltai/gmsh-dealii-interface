// by Raksha Devi

#ifndef GMSH_API_PARALLEL_H
#define GMSH_API_PARALLEL_H

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/point.h>
#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_description.h>

#ifdef DEAL_II_GMSH_WITH_API
#  include <gmsh.h>
#endif

#include <array>
#include <map>
#include <string>
#include <vector>
#include <iostream>

DEAL_II_NAMESPACE_OPEN

namespace GMSH
{
#ifdef DEAL_II_GMSH_WITH_API
template <int dim, int spacedim>
void read_parallel_msh(
    dealii::parallel::fullydistributed::Triangulation<dim, spacedim> &tria,
    const MPI_Comm &mpi_comm,
    const std::string &file_prefix = "mesh_",
    const std::string &file_suffix = ".msh")
{
    using namespace dealii;

    const unsigned int rank = Utilities::MPI::this_mpi_process(mpi_comm);
    const std::string fname = file_prefix + std::to_string(rank) + file_suffix;

    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 0);
    gmsh::open(fname);

    AssertThrow(gmsh::model::getDimension() == dim,
                ExcMessage("Dimension mismatch: Gmsh file dimension " +
                           std::to_string(gmsh::model::getDimension()) +
                           " vs triangulation dimension " + std::to_string(dim)));

    // Retrieve nodes
    std::vector<std::size_t> node_tags;
    std::vector<double> coords, parametric_coords;
    gmsh::model::mesh::getNodes(node_tags, coords, parametric_coords);

    // Prepare triangulation description
    TriangulationDescription::Description<dim, spacedim> triangulation_description;
    triangulation_description.comm = mpi_comm;

    auto &vertices = triangulation_description.coarse_cell_vertices;
    auto &cells = triangulation_description.coarse_cells;
    auto &coarse_cell_ids = triangulation_description.coarse_cell_index_to_coarse_cell_id;
    auto &cell_infos = triangulation_description.cell_infos;


    // Map node tags to indices
    std::map<std::size_t, unsigned int> node_tag_to_index;
    vertices.resize(node_tags.size(), Point<spacedim>());
    for (unsigned int i = 0; i < node_tags.size(); ++i) {
        node_tag_to_index[node_tags[i]] = i;
        for (unsigned int d = 0; d < spacedim; ++d)
            vertices[i][d] = coords[3 * i + d];
    }

    // Retrieve elements
    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities);

    for (const auto &e : entities) {
        const int entity_dim = e.first;
        const int entity_tag = e.second;

        if (entity_dim == dim) {
            std::vector<int> element_types;
            std::vector<std::vector<std::size_t>> element_ids, element_nodes;
            gmsh::model::mesh::getElements(element_types, element_ids, element_nodes, entity_dim, entity_tag);

            for (unsigned int i = 0; i < element_types.size(); ++i) {
                if (element_ids[i].empty()) continue;

                const unsigned int n_vertices = element_nodes[i].size() / element_ids[i].size();

                for (unsigned int j = 0; j < element_ids[i].size(); ++j) {
                    CellData<dim> cell;
                    cell.material_id = 0;

                    for (unsigned int v = 0; v < n_vertices; ++v) {
                        const std::size_t node_tag = element_nodes[i][j * n_vertices + v];
                        cell.vertices[v] = node_tag_to_index[node_tag];
                    }

                    cells.push_back(cell);
                    coarse_cell_ids.push_back(element_ids[i][j]);
                    
                    // // Create cell info with ownership data
                    // TriangulationDescription::Description<dim, spacedim>::CellData info;
                    // info.id = element_ids[i][j];
                    // info.subdomain_id = rank;
                    // info.level_subdomain_id = rank;
                    // cell_infos[0].push_back(info);
                }
            }
        }
    }

    // Use default settings
    triangulation_description.settings = TriangulationDescription::Settings::default_setting;

    // Create the triangulation
    tria.create_triangulation(triangulation_description);

    gmsh::clear();
    gmsh::finalize();
}
#endif
} // namespace GMSH

DEAL_II_NAMESPACE_CLOSE

#endif // GMSH_API_PARALLEL_H
