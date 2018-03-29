/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2016 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, University of Heidelberg, 2000
 */



#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <fstream>
#include <iostream>

// ==========

#include <deal.II/base/mpi.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>

#include <deal.II/distributed/shared_tria.h>

#include <deal.II/grid/filtered_iterator.h>

#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/sparsity_tools.h>

// ==========

using namespace dealii;



template <int dim>
class Step6
{
public:
  Step6 ();
  ~Step6 ();

  void run ();

private:
  void setup_system ();
  void assemble_system ();
  void solve ();
  void refine_grid ();
  void output_results (const unsigned int cycle) const;

  MPI_Comm mpi_communicator;
  const unsigned int n_mpi_processes;
  const unsigned int this_mpi_process;
  ConditionalOStream pcout;

  parallel::shared::Triangulation<dim> triangulation;

  DoFHandler<dim>      dof_handler;
  FE_Q<dim>            fe;
  const QGauss<dim>    quadrature_formula;

  ConstraintMatrix     constraints;

  TrilinosWrappers::SparseMatrix     system_matrix;
  TrilinosWrappers::MPI::Vector      solution;
  TrilinosWrappers::MPI::Vector      system_rhs;

  std::vector<types::global_dof_index> local_dofs_per_process;
  IndexSet                             locally_owned_dofs;
  IndexSet                             locally_relevant_dofs;
};



template <int dim>
double coefficient (const Point<dim> &p)
{
  const Point<dim> centre (0.6,0.6);
  if (p.distance(centre) < 0.2)
    return 5;
  else
    return 1;
}

template <int dim>
double force_function (const Point<dim> &p)
{
  if (std::abs(p[0]*p[1]) > 0.05)
    return 1;
  else
    return -10;
}





template <int dim>
Step6<dim>::Step6 ()
  :
  mpi_communicator (MPI_COMM_WORLD),
  n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
  this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator)),
  pcout (std::cout, this_mpi_process == 0),
  triangulation (mpi_communicator),
  dof_handler (triangulation),
  fe (2),
  quadrature_formula(3)
{}



template <int dim>
Step6<dim>::~Step6 ()
{
  dof_handler.clear ();
}



template <int dim>
void Step6<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  local_dofs_per_process = dof_handler.n_locally_owned_dofs_per_processor();
  locally_owned_dofs = dof_handler.locally_owned_dofs();
  DoFTools::extract_locally_relevant_dofs (dof_handler,locally_relevant_dofs);

  constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
                                           constraints);
  VectorTools::interpolate_boundary_values (dof_handler,
                                            0,
                                            ZeroFunction<dim>(),
                                            constraints);
  constraints.close ();

  DynamicSparsityPattern sparsity_pattern (locally_relevant_dofs);
  DoFTools::make_sparsity_pattern (dof_handler,
                                   sparsity_pattern,
                                   constraints,
                                   /*keep_constrained_dofs =*/ false);
  SparsityTools::distribute_sparsity_pattern (sparsity_pattern,
                                              local_dofs_per_process,
                                              mpi_communicator,
                                              locally_relevant_dofs);
  system_matrix.reinit (locally_owned_dofs,
                        locally_owned_dofs,
                        sparsity_pattern,
                        mpi_communicator);
  solution.reinit (locally_owned_dofs, mpi_communicator);
  system_rhs.reinit (locally_owned_dofs, mpi_communicator);
}



template <int dim>
void Step6<dim>::assemble_system ()
{
  FEValues<dim> fe_values (fe, quadrature_formula,
                           update_values    |  update_gradients |
                           update_quadrature_points  |  update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       cell_rhs (dofs_per_cell);

  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
  cell (IteratorFilters::LocallyOwnedCell(),
      dof_handler.begin_active()),
  endc (IteratorFilters::LocallyOwnedCell(),
      dof_handler.end());
  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      cell_rhs = 0;

      fe_values.reinit (cell);

      for (unsigned int q_index=0; q_index<n_q_points; ++q_index)
        {
          const double current_coefficient = coefficient<dim>
                                             (fe_values.quadrature_point (q_index));
          const double current_force_function = force_function<dim>
                                                (fe_values.quadrature_point (q_index));
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              for (unsigned int j=0; j<dofs_per_cell; ++j)
                cell_matrix(i,j) += (current_coefficient *
                                     fe_values.shape_grad(i,q_index) *
                                     fe_values.shape_grad(j,q_index) *
                                     fe_values.JxW(q_index));

              cell_rhs(i) += (fe_values.shape_value(i,q_index) *
                              current_force_function *
                              fe_values.JxW(q_index));
            }
        }

      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global (cell_matrix,
                                              cell_rhs,
                                              local_dof_indices,
                                              system_matrix,
                                              system_rhs);
    }

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
}




template <int dim>
void Step6<dim>::solve ()
{
  SolverControl                            solver_control (1000, 1e-12);
  SolverCG<TrilinosWrappers::MPI::Vector>  solver (solver_control);

  TrilinosWrappers::PreconditionSSOR preconditioner;
  preconditioner.initialize(system_matrix);

  solver.solve (system_matrix, solution, system_rhs,
                preconditioner);

  constraints.distribute (solution);
}



template <int dim>
void Step6<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  TrilinosWrappers::MPI::Vector locally_relevant_solution;
  locally_relevant_solution.reinit (locally_owned_dofs,
                                    locally_relevant_dofs,
                                    mpi_communicator);
  locally_relevant_solution = solution;

  KellyErrorEstimator<dim>::estimate (dof_handler,
                                      QGauss<dim-1>(3),
                                      typename FunctionMap<dim>::type(),
                                      locally_relevant_solution,
                                      estimated_error_per_cell);

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                   estimated_error_per_cell,
                                                   0.3, 0.03);

  triangulation.execute_coarsening_and_refinement ();
}



template <int dim>
void Step6<dim>::output_results (const unsigned int cycle) const
{
  Assert (cycle < 10, ExcNotImplemented());

  std::string filename = "grid-";
  filename += ('0' + cycle);
  filename += ".vtu";

  DataOut<dim> data_out;

  TrilinosWrappers::MPI::Vector locally_relevant_solution;
  locally_relevant_solution.reinit (locally_owned_dofs,
                                    locally_relevant_dofs,
                                    mpi_communicator);
  locally_relevant_solution = solution;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (locally_relevant_solution, "solution");

  Vector<double> subdomain_id (triangulation.n_active_cells());
  typename Triangulation<dim>::active_cell_iterator
  cell = triangulation.begin_active(),
  endc = triangulation.end();
  unsigned int cell_no=0;
  for (; cell!=endc; ++cell, ++cell_no)
    {
      subdomain_id[cell_no] = cell->subdomain_id();
    }
  data_out.add_data_vector (subdomain_id, "subdomain_id");

  data_out.build_patches (fe.tensor_degree());
  data_out.write_vtu_in_parallel (filename.c_str(), mpi_communicator);
}



template <int dim>
void Step6<dim>::run ()
{
  for (unsigned int cycle=0; cycle<8; ++cycle)
    {
      pcout << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
        {
          GridGenerator::hyper_cube (triangulation);
          triangulation.refine_global (3);
        }
      else
        refine_grid ();


      pcout << "   Number of active cells:       "
                << triangulation.n_active_cells()
                << " (by partition:";
      for (unsigned int p=0; p<n_mpi_processes; ++p)
        pcout << (p==0 ? ' ' : '+')
              << (GridTools::count_cells_with_subdomain_association (triangulation,p));
      pcout << ")" << std::endl;

      setup_system ();

      pcout << "   Number of degrees of freedom: "
                << dof_handler.n_dofs()
                << " (by partition:";
      for (unsigned int p=0; p<n_mpi_processes; ++p)
        pcout << (p==0 ? ' ' : '+')
              << (DoFTools::
                  count_dofs_with_subdomain_association (dof_handler,p));
      pcout << ")" << std::endl;

      assemble_system ();
      solve ();
      output_results (cycle);
    }
}



int main (int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  try
    {
      Step6<2> laplace_problem_2d;
      laplace_problem_2d.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
