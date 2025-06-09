#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>
#include <deal.II/base/mpi.h>

#include <fstream>
#include "gmsh_api_parallel.h"
#include "tests.h"

auto main(int argc, char **argv) -> int
{
  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm mpi_comm = MPI_COMM_WORLD;

  // Get the rank and number of processes
  const unsigned int rank = dealii::Utilities::MPI::this_mpi_process(mpi_comm);
  const unsigned int n_processes = dealii::Utilities::MPI::n_mpi_processes(mpi_comm);

  // Ensure we have exactly 2 processes
  if (n_processes != 2)
  {
    if (rank == 0)
      std::cerr << "This program requires exactly 2 MPI processes." << std::endl;
    MPI_Finalize();
    return 1;
  }

  // Initialize logging
  initlog(rank);

  // Create a 2D triangulation
  dealii::Triangulation<2, 2> triangulation;

  // Define mesh files for each process
  std::string mesh_files[2] = {
      SOURCE_DIR "/../grids/new_meshes/test_grid_2_1.msh",
      SOURCE_DIR "/../grids/new_meshes/test_grid_2_2.msh"};

  // Each process reads its corresponding mesh file
  GMSH::read_msh(triangulation, mesh_files[rank]);

  // Output the number of cells and vertices for this process
  deallog << "Process " << rank << ": Number of cells = " << triangulation.n_active_cells() << std::endl;
  deallog << "Process " << rank << ": Number of vertices = " << triangulation.n_vertices() << std::endl;

  // Finalize MPI
  MPI_Finalize();
  return 0;
}
