#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <ostream>

// #include "gmsh_api_parallel.h"
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

  // Create the fully distributed triangulation
  parallel::fullydistributed::Triangulation<2, 2> tria(mpi_comm);

  parallel::distributed::Triangulation<2, 2> distributed_tria(mpi_comm);
  GridGenerator::hyper_cube(distributed_tria, -1, 1);
  distributed_tria.refine_global(3);

  auto description =
    TriangulationDescription::Utilities::create_description_from_triangulation(
      distributed_tria, mpi_comm);


  tria.create_triangulation(description);
  // Ordered output

  deallog << "Coarse cells: " << description.coarse_cells.size() << std::endl
          << "Coarse vertices: " << description.coarse_cell_vertices.size()
          << std::endl
          << "Coarse cell index to coarse cell id: "
          << description.coarse_cell_index_to_coarse_cell_id.size() << std::endl
          << "Cell infos: " << description.cell_infos.size() << std::endl;

  deallog << "Rank " << rank << " owns " << tria.n_active_cells()
          << " cells and " << tria.n_vertices() << " vertices." << std::endl;



  return 0;
}