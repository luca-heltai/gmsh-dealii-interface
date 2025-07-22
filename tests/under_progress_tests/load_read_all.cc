#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>

#include <fstream>
#include <string>
#include <vector>

#include "gmsh_api.h"
#include "tests.h"

auto
main() -> int
{
  initlog(1);

  // List of split mesh files
  std::vector<std::string> split_mesh_files = {
    SOURCE_DIR "/../grids/new_meshes/test_grid_2_1.msh",
    SOURCE_DIR "/../grids/new_meshes/test_grid_2_2.msh",
    SOURCE_DIR "/../grids/new_meshes/test_grid_2_3.msh"};

  // Read and output details for each split mesh file
  for (const auto &mesh_file : split_mesh_files)
    {
      // Create a 2D triangulation for each mesh file
      dealii::Triangulation<2, 2> triangulation;

      // Open and read the mesh file
      GMSH::read_msh(triangulation, mesh_file);

      // Output the number of cells and vertices for the current mesh
      deallog << "File: " << mesh_file << std::endl;
      deallog << "Number of cells: " << triangulation.n_active_cells()
              << std::endl;
      deallog << "Number of vertices: " << triangulation.n_vertices()
              << std::endl;
    }

  // Initialize Gmsh API for merging
  gmsh::initialize();

  // Merge all split mesh files
  for (const auto &file : split_mesh_files)
    {
      std::cout << "Merging file: " << file << std::endl;
      gmsh::merge(file);
    }

  // Optionally remove duplicate nodes
  gmsh::model::mesh::removeDuplicateNodes();

  // Save the merged mesh to a new file
  std::string merged_mesh_file =
    SOURCE_DIR "/../grids/new_meshes/test_grid_2_merged.msh";
  gmsh::write(merged_mesh_file);
  std::cout << "Merged mesh saved to: " << merged_mesh_file << std::endl;

  // Finalize Gmsh API
  gmsh::finalize();

  // Read and output details for the merged mesh file
  dealii::Triangulation<2, 2> merged_triangulation;

  // Open and read the merged mesh file
  GMSH::read_msh(merged_triangulation, merged_mesh_file);

  // Output the number of cells and vertices in the merged mesh
  deallog << "Merged File: " << merged_mesh_file << std::endl;
  deallog << "Number of cells in merged mesh: "
          << merged_triangulation.n_active_cells() << std::endl;
  deallog << "Number of vertices in merged mesh: "
          << merged_triangulation.n_vertices() << std::endl;

  return 0;
}