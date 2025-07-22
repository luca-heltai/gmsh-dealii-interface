//

#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/grid/grid_out.h>

#include <fstream>

#include "gmsh_api_parallel.h"
#include "tests.h"
using namespace dealii;



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    mpilog;
  MPI_Comm                         mpi_comm = MPI_COMM_WORLD;

  int rank, size;
  MPI_Comm_rank(mpi_comm, &rank);
  MPI_Comm_size(mpi_comm, &size);


  parallel::fullydistributed::Triangulation<2, 2> tria(mpi_comm);


  GMSH::read_partitioned_msh(tria,
                             mpi_comm,
                             SOURCE_DIR "/../grids/add_meshes/unit-square");



  deallog << "Rank " << rank << " owns " << tria.n_active_cells()
          << " cells and " << tria.n_vertices() << " vertices." << std::endl;

  unsigned int n_locally_owned_cells = 0;
  unsigned int n_ghost_cells         = 0;

  for (const auto &cell : tria.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        ++n_locally_owned_cells;
      else if (cell->is_ghost())
        ++n_ghost_cells;
    }

  deallog << "Rank " << rank << " owns " << n_locally_owned_cells
          << " cells and has " << n_ghost_cells << " ghost cells." << std::endl;

  std::cout << "Rank " << rank << " finished deallog output" << std::endl;

  std::cout << "Rank " << rank << " exiting main()" << std::endl;

  return 0;
}
