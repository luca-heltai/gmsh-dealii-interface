/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2021 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// Step 40 on partitioned quad meshes.
// Modifications with respect to original step-40:
// - grid generation via GridGenerator::subdivided_hyper_rectangle and
//   GridGenerator::subdivided_hyper_rectangle_with_quads (instead of
//   GridGenerator::hyper_cube);
// - removed local, adaptive refinement.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/generic_linear_algebra.h>

namespace LA
{
#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \
  !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
  using namespace dealii::LinearAlgebraPETSc;
#  define USE_PETSC_LA
#elif defined(DEAL_II_WITH_TRILINOS)
  using namespace dealii::LinearAlgebraTrilinos;
#else
#  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
#endif
} // namespace LA

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/fully_distributed_tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "gmsh_api_parallel.h"
#include "tests.h"

namespace Step40
{
  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem();
    void
    run();

  private:
    void
    setup_system();
    void
    assemble_system();
    void
    solve();
    void
    output_results() const;

    MPI_Comm                                       mpi_communicator;
    MappingQ<dim>                                  mapping;
    parallel::fullydistributed::Triangulation<dim> triangulation;
    FE_Q<dim>                                      fe;
    DoFHandler<dim>                                dof_handler;

    IndexSet locally_owned_dofs;
    IndexSet locally_relevant_dofs;

    AffineConstraints<double> constraints;
    LA::MPI::SparseMatrix     system_matrix;
    LA::MPI::Vector           locally_relevant_solution;
    LA::MPI::Vector           system_rhs;
  };

  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem()
    : mpi_communicator(MPI_COMM_WORLD)
    , mapping(1) // linear mapping for FE_Q(2) can also be 1 here
    , triangulation(mpi_communicator)
    , fe(1)
    , dof_handler(triangulation)
  {}

  template <int dim>
  void
  LaplaceProblem<dim>::setup_system()
  {
    std::cout << "Rank " << Utilities::MPI::this_mpi_process(mpi_communicator)
              << ": Starting setup_system()" << std::endl;

    dof_handler.distribute_dofs(fe);
    std::cout << "Rank " << Utilities::MPI::this_mpi_process(mpi_communicator)
              << ": DOFs distributed successfully" << std::endl;

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    std::cout << "Rank " << Utilities::MPI::this_mpi_process(mpi_communicator)
              << ": Locally owned DOFs: " << locally_owned_dofs.n_elements()
              << std::endl;

    locally_relevant_dofs =
      DoFTools::extract_locally_relevant_dofs(dof_handler);
    std::cout << "Rank " << Utilities::MPI::this_mpi_process(mpi_communicator)
              << ": Locally relevant DOFs: "
              << locally_relevant_dofs.n_elements() << std::endl;

    MPI_Barrier(mpi_communicator);

    locally_relevant_solution.reinit(locally_owned_dofs,
                                     locally_relevant_dofs,
                                     mpi_communicator);
    std::cout << "Rank " << Utilities::MPI::this_mpi_process(mpi_communicator)
              << ": Solution vector reinitialized" << std::endl;

    system_rhs.reinit(locally_owned_dofs, mpi_communicator);
    std::cout << "Rank " << Utilities::MPI::this_mpi_process(mpi_communicator)
              << ": RHS vector reinitialized" << std::endl;

    constraints.clear();
    constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
    std::cout << "Rank " << Utilities::MPI::this_mpi_process(mpi_communicator)
              << ": Constraints reinitialized" << std::endl;
    DoFTools::make_hanging_node_constraints(dof_handler, constraints);
    VectorTools::interpolate_boundary_values(
      mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints);
    constraints.close();

    DynamicSparsityPattern dsp(locally_relevant_dofs);
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
    SparsityTools::distribute_sparsity_pattern(dsp,
                                               dof_handler.locally_owned_dofs(),
                                               mpi_communicator,
                                               locally_relevant_dofs);

    system_matrix.reinit(locally_owned_dofs,
                         locally_owned_dofs,
                         dsp,
                         mpi_communicator);
  }

  template <int dim>
  void
  LaplaceProblem<dim>::assemble_system()
  {
    const QGauss<dim> quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(mapping,
                            fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators() |
                              IteratorFilters::LocallyOwnedCell())
      {
        cell_matrix = 0.;
        cell_rhs    = 0.;
        fe_values.reinit(cell);

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            const double rhs_value =
              std::sin(numbers::PI * fe_values.quadrature_point(q_point)[0]);

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  cell_matrix(i, j) += fe_values.shape_grad(i, q_point) *
                                       fe_values.shape_grad(j, q_point) *
                                       fe_values.JxW(q_point);

                cell_rhs(i) += rhs_value * fe_values.shape_value(i, q_point) *
                               fe_values.JxW(q_point);
              }
          }

        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
      }

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }

  template <int dim>
  void
  LaplaceProblem<dim>::solve()
  {
    LA::MPI::Vector completely_distributed_solution(locally_owned_dofs,
                                                    mpi_communicator);
    SolverControl   solver_control(dof_handler.n_dofs(), 1e-10);

#ifdef USE_PETSC_LA
    LA::SolverCG solver(solver_control, mpi_communicator);
#else
    LA::SolverCG solver(solver_control);
#endif

    LA::MPI::PreconditionAMG                 preconditioner;
    LA::MPI::PreconditionAMG::AdditionalData data;

#ifdef USE_PETSC_LA
    data.symmetric_operator = true;
#endif

    preconditioner.initialize(system_matrix, data);

    solver.solve(system_matrix,
                 completely_distributed_solution,
                 system_rhs,
                 preconditioner);

    constraints.distribute(completely_distributed_solution);

    locally_relevant_solution = completely_distributed_solution;
  }

  template <int dim>
  void
  LaplaceProblem<dim>::output_results() const
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(locally_relevant_solution, "solution");
    data_out.build_patches(mapping, fe.degree);

    const std::string filename =
      "solution-" +
      Utilities::int_to_string(
        Utilities::MPI::this_mpi_process(mpi_communicator), 4) +
      ".vtu";
    std::ofstream output(filename);
    data_out.write_vtu(output);
  }

  template <int dim>
  void
  LaplaceProblem<dim>::run()
  {
    deallog << "Running on "
            << Utilities::MPI::n_mpi_processes(mpi_communicator)
            << " MPI rank(s)..." << std::endl;
    // GMSH::read_partitioned_msh(triangulation, mpi_communicator,
    //                            "/home/raksha/dealii/parallel_gmsh_codes/grids/add_meshes/unit-square");
    GMSH::read_partitioned_msh(triangulation,
                               mpi_communicator,
                               SOURCE_DIR "/../grids/square-quad",
                               ".msh");

    deallog << "Triangulation has " << triangulation.n_active_cells()
            << " active cells." << std::endl;
    std::cout << "Rank " << Utilities::MPI::this_mpi_process(mpi_communicator)
              << ": Local cells = " << triangulation.n_active_cells()
              << ", Global cells = " << triangulation.n_global_active_cells()
              << std::endl;

    // Count and print cells owned and ghosts
    unsigned int ghost_cells = 0, owned_cells = 0;
    for (const auto &cell : triangulation.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          ++owned_cells;
        else if (cell->is_ghost())
          ++ghost_cells;
      }
    std::cout << "Rank " << Utilities::MPI::this_mpi_process(mpi_communicator)
              << ": Owned cells: " << owned_cells
              << ", Ghost cells: " << ghost_cells << std::endl;

    setup_system();
    assemble_system();
    solve();
    output_results();
  }

} // namespace Step40

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      Step40::LaplaceProblem<2>        laplace_problem;
      laplace_problem.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl;
      std::cerr << "Aborting!" << std::endl;
      std::cerr << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl;
      std::cerr << "Aborting!" << std::endl;
      std::cerr << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;
}
