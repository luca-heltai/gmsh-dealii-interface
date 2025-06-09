
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>

#include <fstream>

#include "gmsh_api.h"
#include "tests.h"

auto
main() -> int
{
  initlog(1);

  // Create a 2D triangulation
  dealii::Triangulation<2, 2> triangulation;

  // Open and read the mesh file
  std::string mesh_file = SOURCE_DIR "/../grids/new_meshes/test_grid_2.msh";

  GMSH::read_msh(triangulation, mesh_file);

  // Output the number of cells and vertices
  deallog << "Number of cells: " << triangulation.n_active_cells() << std::endl;
  deallog << "Number of vertices: " << triangulation.n_vertices() << std::endl;

  return 0;
}
