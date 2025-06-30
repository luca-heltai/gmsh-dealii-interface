#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>
#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/grid/grid_out.h>
#include "gmsh_api_parallel.h"
#include "tests.h"
using namespace dealii;

int main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll mpilog;
  MPI_Comm mpi_comm = MPI_COMM_WORLD;

  int rank, size;
  MPI_Comm_rank(mpi_comm, &rank);
  MPI_Comm_size(mpi_comm, &size);

  // Create the fully distributed triangulation
  parallel::fullydistributed::Triangulation<2, 2> tria(mpi_comm);

  // Read the partitioned mesh 
  GMSH::read_parallel_msh(tria,
                          mpi_comm,
                          SOURCE_DIR "/../grids/add_meshes/unit-square_");

                          
  // Ordered output


   deallog << "Rank " << rank << " owns " << tria.n_active_cells()
               << " cells and " << tria.n_vertices() << " vertices." << std::endl;



  return 0;
}