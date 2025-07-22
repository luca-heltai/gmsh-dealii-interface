//  by Raksha Devi

#include <deal.II/base/point.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>

#include <gmsh.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "gmsh_api.h"
#include "tests.h"
int
main(int argc, char **argv)
{
  // Initialize the Gmsh API
  gmsh::initialize();

  // List of partitioned mesh files with directory path
  std::vector<std::string> partitioned_files = {
    SOURCE_DIR "/../grids/new_meshes/test_grid_2_1.msh",
    SOURCE_DIR "/../grids/new_meshes/test_grid_2_2.msh",
    SOURCE_DIR "/../grids/new_meshes/test_grid_2_3.msh",
  };

  // Merge all partitioned meshes
  for (const auto &file : partitioned_files)
    {
      gmsh::merge(file);
    }

  // remove duplicate nodes and physical  groups
  gmsh::model::mesh::removeDuplicateNodes();
  gmsh::model::removePhysicalGroups();
  // Create a 2D triangulation
  // dealii::Triangulation<2, 2> triangulation;

  // Save the merged mesh to a new file
  std::string output_file =
    SOURCE_DIR "/../grids/new_meshes/test_grid_2_merged.msh";
  gmsh::write(output_file);

  // Get nodes and elements
  std::vector<std::size_t> node_tags;
  std::vector<double>      node_coords;
  std::vector<double>      parametric_coords;

  gmsh::model::mesh::getNodes(
    node_tags, node_coords, parametric_coords, -1, -1, true, false);

  // Output the number of vertices
  std::cout << "DEAL::Number of vertices: " << node_tags.size() << std::endl;

  // Finalize Gmsh API
  gmsh::finalize();


  // Count the number of cells and vertices in the merged mesh
  dealii::Triangulation<2, 2> triangulation;

  // Read the merged mesh file
  dealii::GridIn<2> grid_in;
  grid_in.attach_triangulation(triangulation);
  std::ifstream input_file(SOURCE_DIR
                           "/../grids/new_meshes/test_grid_2_merged.msh");
  grid_in.read_msh(input_file);

  // Get nodes and elements    std::vector<double> node_coords;
  // gmsh::model::mesh::getNodes(node_tags, node_coords);
  // int num_vertices = node_tags.size();


  // Output the number of cells
  std::cout << "DEAL::Number of cells: " << triangulation.n_active_cells()
            << std::endl;

  return 0;
}