// by Raksha Devi

#include <deal.II/base/mpi.h>
#include <deal.II/distributed/fully_distributed_tria.h>
#include "gmsh_api_parallel.h"

using namespace dealii;

int main(int argc, char **argv)
{
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    MPI_Comm mpi_comm = MPI_COMM_WORLD;

    parallel::fullydistributed::Triangulation<2, 2> tria(mpi_comm);

    GMSH::merge_meshes_and_create_triangulation(tria,
                                                mpi_comm,
                                                SOURCE_DIR "/../grids/New_Meshes/test_grid_2_",
                                                ".msh");

    return 0;
}