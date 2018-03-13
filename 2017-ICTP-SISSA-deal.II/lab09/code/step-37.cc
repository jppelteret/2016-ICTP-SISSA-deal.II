/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2017 by the deal.II authors
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
 * Authors: Martin Kronbichler, based on step-37
 */


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>

#include <deal.II/numerics/vector_tools.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <iostream>
#include <fstream>
#include <sstream>


namespace Step37
{
  using namespace dealii;


  const unsigned int degree_finite_element = 2;
  const unsigned int dimension = 3;



  template <int dim, int fe_degree, typename number>
  class LaplaceOperator : public MatrixFreeOperators::Base<dim,LinearAlgebra::distributed::Vector<number> >
  {
  public:
    typedef number value_type;

    LaplaceOperator ();

    virtual void compute_diagonal();

  private:
    virtual void apply_add(LinearAlgebra::distributed::Vector<number> &dst,
                           const LinearAlgebra::distributed::Vector<number> &src) const;

    void local_apply (const MatrixFree<dim,number>                     &data,
                      LinearAlgebra::distributed::Vector<number>       &dst,
                      const LinearAlgebra::distributed::Vector<number> &src,
                      const std::pair<unsigned int,unsigned int>       &cell_range) const;

    void local_compute_diagonal (const MatrixFree<dim,number>                     &data,
                                 LinearAlgebra::distributed::Vector<number>       &dst,
                                 const unsigned int                               &dummy,
                                 const std::pair<unsigned int,unsigned int>       &cell_range) const;
  };



  template <int dim, int fe_degree, typename number>
  LaplaceOperator<dim,fe_degree,number>::LaplaceOperator ()
    :
    MatrixFreeOperators::Base<dim, LinearAlgebra::distributed::Vector<number> >()
  {}




  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>
  ::local_apply (const MatrixFree<dim,number>                     &data,
                 LinearAlgebra::distributed::Vector<number>       &dst,
                 const LinearAlgebra::distributed::Vector<number> &src,
                 const std::pair<unsigned int,unsigned int>       &cell_range) const
  {
    FEEvaluation<dim,fe_degree,fe_degree+1,1,number> phi (data);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        phi.reinit (cell);
        phi.read_dof_values(src);
        phi.evaluate (false, true);
        for (unsigned int q=0; q<phi.n_q_points; ++q)
          phi.submit_gradient (phi.get_gradient(q), q);
        phi.integrate (false, true);
        phi.distribute_local_to_global (dst);
      }
  }



  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>
  ::apply_add (LinearAlgebra::distributed::Vector<number>       &dst,
               const LinearAlgebra::distributed::Vector<number> &src) const
  {
    this->data->cell_loop (&LaplaceOperator::local_apply, this, dst, src);
  }



  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>
  ::compute_diagonal ()
  {
    this->inverse_diagonal_entries.
    reset(new DiagonalMatrix<LinearAlgebra::distributed::Vector<number> >());
    LinearAlgebra::distributed::Vector<number> &inverse_diagonal =
      this->inverse_diagonal_entries->get_vector();
    this->data->initialize_dof_vector(inverse_diagonal);
    unsigned int dummy = 0;
    this->data->cell_loop (&LaplaceOperator::local_compute_diagonal, this,
                           inverse_diagonal, dummy);

    this->set_constrained_entries_to_one(inverse_diagonal);

    for (unsigned int i=0; i<inverse_diagonal.local_size(); ++i)
      {
        Assert(inverse_diagonal.local_element(i) > 0.,
               ExcMessage("No diagonal entry in a positive definite operator "
                          "should be zero"));
        inverse_diagonal.local_element(i) =
          1./inverse_diagonal.local_element(i);
      }
  }



  template <int dim, int fe_degree, typename number>
  void
  LaplaceOperator<dim,fe_degree,number>
  ::local_compute_diagonal (const MatrixFree<dim,number>               &data,
                            LinearAlgebra::distributed::Vector<number> &dst,
                            const unsigned int &,
                            const std::pair<unsigned int,unsigned int> &cell_range) const
  {
    FEEvaluation<dim,fe_degree,fe_degree+1,1,number> phi (data);

    AlignedVector<VectorizedArray<number> > diagonal(phi.dofs_per_cell);

    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell)
      {
        phi.reinit (cell);
        for (unsigned int i=0; i<phi.dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<phi.dofs_per_cell; ++j)
              phi.submit_dof_value(VectorizedArray<number>(), j);
            phi.submit_dof_value(make_vectorized_array<number>(1.), i);

            phi.evaluate (false, true);
            for (unsigned int q=0; q<phi.n_q_points; ++q)
              phi.submit_gradient (phi.get_gradient(q), q);
            phi.integrate (false, true);
            diagonal[i] = phi.get_dof_value(i);
          }
        for (unsigned int i=0; i<phi.dofs_per_cell; ++i)
          phi.submit_dof_value(diagonal[i], i);
        phi.distribute_local_to_global (dst);
      }
  }




  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem ();
    void run ();

  private:
    void setup_system ();
    void assemble_rhs ();
    void solve ();
    void compare_performance_matrix_vector () const;

#ifdef DEAL_II_WITH_P4EST
    parallel::distributed::Triangulation<dim>  triangulation;
#else
    Triangulation<dim>                         triangulation;
#endif

    FE_Q<dim>                                  fe;
    DoFHandler<dim>                            dof_handler;

    ConstraintMatrix                           constraints;
    typedef LaplaceOperator<dim,degree_finite_element,double> SystemMatrixType;
    SystemMatrixType                           system_matrix;

    MGConstrainedDoFs                          mg_constrained_dofs;
    typedef LaplaceOperator<dim,degree_finite_element,float>  LevelMatrixType;
    MGLevelObject<LevelMatrixType>             mg_matrices;

    LinearAlgebra::distributed::Vector<double> solution;
    LinearAlgebra::distributed::Vector<double> system_rhs;

    ConditionalOStream                         pcout;
    mutable TimerOutput                        timer_output;
  };



  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem ()
    :
#ifdef DEAL_II_WITH_P4EST
    triangulation(MPI_COMM_WORLD,
                  Triangulation<dim>::limit_level_difference_at_vertices,
                  parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy),
#else
    triangulation (Triangulation<dim>::limit_level_difference_at_vertices),
#endif
    fe (degree_finite_element),
    dof_handler (triangulation),
    pcout (std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0),
    timer_output (pcout, TimerOutput::summary, TimerOutput::wall_times)
  {}





  template <int dim>
  void LaplaceProblem<dim>::setup_system ()
  {
    {
      TimerOutput::Scope t(timer_output, "Distribute DoFs & B.C.");

      system_matrix.clear();
      mg_matrices.clear_elements();

      dof_handler.distribute_dofs (fe);
      dof_handler.distribute_mg_dofs (fe);

      IndexSet locally_relevant_dofs;
      DoFTools::extract_locally_relevant_dofs (dof_handler,
                                               locally_relevant_dofs);

      constraints.clear();
      constraints.reinit(locally_relevant_dofs);
      DoFTools::make_hanging_node_constraints(dof_handler,
                                              constraints);
      VectorTools::interpolate_boundary_values (dof_handler,
                                                0,
                                                ZeroFunction<dim>(),
                                                constraints);
      constraints.close();
    }

    {
      TimerOutput::Scope t(timer_output, "Setup matrix-free system");
      typename MatrixFree<dim,double>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
        MatrixFree<dim,double>::AdditionalData::none;
      additional_data.mapping_update_flags = (update_gradients | update_JxW_values |
                                              update_quadrature_points);
      std::shared_ptr<MatrixFree<dim,double> >
      system_mf_storage(new MatrixFree<dim,double>());
      system_mf_storage->reinit (dof_handler, constraints, QGauss<1>(fe.degree+1),
                                 additional_data);
      system_matrix.initialize (system_mf_storage);

      system_matrix.initialize_dof_vector(solution);
      system_matrix.initialize_dof_vector(system_rhs);
    }

    TimerOutput::Scope t(timer_output, "Setup matrix-free levels");

    const unsigned int nlevels = triangulation.n_global_levels();
    mg_matrices.resize(0, nlevels-1);

    std::set<types::boundary_id> dirichlet_boundary;
    dirichlet_boundary.insert(0);
    mg_constrained_dofs.initialize(dof_handler);
    mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, dirichlet_boundary);

    for (unsigned int level=0; level<nlevels; ++level)
      {
        IndexSet relevant_dofs;
        DoFTools::extract_locally_relevant_level_dofs(dof_handler, level,
                                                      relevant_dofs);
        ConstraintMatrix level_constraints;
        level_constraints.reinit(relevant_dofs);
        level_constraints.add_lines(mg_constrained_dofs.get_boundary_indices(level));
        level_constraints.close();

        typename MatrixFree<dim,float>::AdditionalData additional_data;
        additional_data.tasks_parallel_scheme =
          MatrixFree<dim,float>::AdditionalData::none;
        additional_data.mapping_update_flags = (update_gradients | update_JxW_values |
                                                update_quadrature_points);
        additional_data.level_mg_handler = level;
        std::shared_ptr<MatrixFree<dim,float> >
        mg_mf_storage_level(new MatrixFree<dim,float>());
        mg_mf_storage_level->reinit(dof_handler, level_constraints,
                                    QGauss<1>(fe.degree+1), additional_data);

        mg_matrices[level].initialize(mg_mf_storage_level, mg_constrained_dofs,
                                      level);
      }
  }




  template <int dim>
  void LaplaceProblem<dim>::assemble_rhs ()
  {
    TimerOutput::Scope t(timer_output, "Assemble rhs");
    system_rhs = 0;
    FEEvaluation<dim,degree_finite_element> phi(*system_matrix.get_matrix_free());
    for (unsigned int cell=0; cell<system_matrix.get_matrix_free()->n_macro_cells();
         ++cell)
      {
        phi.reinit(cell);
        for (unsigned int q=0; q<phi.n_q_points; ++q)
          phi.submit_value(make_vectorized_array<double>(1.0), q);
        phi.integrate(true, false);
        phi.distribute_local_to_global(system_rhs);
      }
    system_rhs.compress(VectorOperation::add);
  }




  template <int dim>
  void LaplaceProblem<dim>::solve ()
  {
    MGTransferMatrixFree<dim,float> mg_transfer(mg_constrained_dofs);
    {
      TimerOutput::Scope t(timer_output, "MG build transfer");
      mg_transfer.build(dof_handler);
    }

    typedef PreconditionChebyshev<LevelMatrixType,LinearAlgebra::distributed::Vector<float> > SmootherType;
    mg::SmootherRelaxation<SmootherType, LinearAlgebra::distributed::Vector<float> >
      mg_smoother;
    {
      TimerOutput::Scope t(timer_output, "MG build smoother");

      MGLevelObject<typename SmootherType::AdditionalData> smoother_data;
      smoother_data.resize(0, triangulation.n_global_levels()-1);
      for (unsigned int level = 0; level<triangulation.n_global_levels(); ++level)
        {
          if (level > 0)
            {
              smoother_data[level].smoothing_range = 15.;
              smoother_data[level].degree = 4;
              smoother_data[level].eig_cg_n_iterations = 10;
            }
          else
            {
              smoother_data[0].smoothing_range = 1e-3;
              smoother_data[0].degree = numbers::invalid_unsigned_int;
              smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m();
            }
          mg_matrices[level].compute_diagonal();
          smoother_data[level].preconditioner = mg_matrices[level].get_matrix_diagonal_inverse();
        }
      mg_smoother.initialize(mg_matrices, smoother_data);
    }

    MGCoarseGridApplySmoother<LinearAlgebra::distributed::Vector<float> > mg_coarse;
    mg_coarse.initialize(mg_smoother);

    mg::Matrix<LinearAlgebra::distributed::Vector<float> > mg_matrix(mg_matrices);

    MGLevelObject<MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType> > mg_interface_matrices;
    mg_interface_matrices.resize(0, triangulation.n_global_levels()-1);
    for (unsigned int level=0; level<triangulation.n_global_levels(); ++level)
      mg_interface_matrices[level].initialize(mg_matrices[level]);
    mg::Matrix<LinearAlgebra::distributed::Vector<float> > mg_interface(mg_interface_matrices);

    Multigrid<LinearAlgebra::distributed::Vector<float> > mg(mg_matrix,
                                                             mg_coarse,
                                                             mg_transfer,
                                                             mg_smoother,
                                                             mg_smoother);
    mg.set_edge_matrices(mg_interface, mg_interface);

    PreconditionMG<dim, LinearAlgebra::distributed::Vector<float>,
                   MGTransferMatrixFree<dim,float> >
                   preconditioner(dof_handler, mg, mg_transfer);


    SolverControl solver_control (100, 1e-12*system_rhs.l2_norm());
    SolverCG<LinearAlgebra::distributed::Vector<double> > cg (solver_control);

    TimerOutput::Scope t(timer_output, "Solve system");
    Timer time;
    cg.solve (system_matrix, solution, system_rhs,
              preconditioner);

    constraints.distribute(solution);

    pcout << "Solve for " << std::setw(8) << dof_handler.n_dofs()
          << " DoFs in "
          << solver_control.last_step()
          << " CG iterations, solver time = "
          << time.wall_time() << "s" << std::endl;
  }



  template <int dim>
  void LaplaceProblem<dim>::compare_performance_matrix_vector() const
  {
    if ((double)dof_handler.n_dofs() * fe.dofs_per_cell
        > 2e8)
      {
        pcout << "Sparse matrix would take more than approx. 3 GB, skip test"
              << std::endl;
        return;
      }

    TrilinosWrappers::SparseMatrix sparse_matrix;
    IndexSet locally_relevant_dofs = constraints.get_local_lines();
    {
      TimerOutput::Scope t(timer_output, "Create sparse matrix");
      DynamicSparsityPattern dsp(locally_relevant_dofs);
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
      SparsityTools::distribute_sparsity_pattern(dsp,
                                                 dof_handler.n_locally_owned_dofs_per_processor(),
                                                 MPI_COMM_WORLD,
                                                 locally_relevant_dofs);
      sparse_matrix.reinit(dof_handler.locally_owned_dofs(),
                           dof_handler.locally_owned_dofs(),
                           dsp, MPI_COMM_WORLD);
    }
    {
      TimerOutput::Scope t(timer_output, "Assemble sparse matrix");

      const QGauss<dim>  quadrature_formula(fe.degree+1);
      FEValues<dim> fe_values (fe, quadrature_formula,
                               update_values    |  update_gradients |
                               update_JxW_values);

      const unsigned int   dofs_per_cell = fe.dofs_per_cell;
      const unsigned int   n_q_points    = quadrature_formula.size();
      FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
      std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

      typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
      for (; cell!=endc; ++cell)
        if (cell->is_locally_owned())
          {
            cell_matrix = 0;

            fe_values.reinit (cell);
            for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                  cell_matrix(i,j) += (fe_values.shape_grad(i,q_point) *
                                       fe_values.shape_grad(j,q_point) *
                                       fe_values.JxW(q_point));
            cell->get_dof_indices (local_dof_indices);
            constraints.distribute_local_to_global (cell_matrix,
                                                    local_dof_indices,
                                                    sparse_matrix);
          }
      sparse_matrix.compress (VectorOperation::add);
    }
    LinearAlgebra::distributed::Vector<double> v1, v2;
    v1.reinit(solution);
    v2.reinit(solution);

    const unsigned int n_products = 50;
    Timer time;
    for (unsigned int i=0; i<n_products; ++i)
      system_matrix.vmult(v1, solution);
    const double time_mf = time.wall_time()/n_products;

    time.restart();
    for (unsigned int i=0; i<n_products; ++i)
      sparse_matrix.vmult(v2, solution);
    const double time_sparse = time.wall_time()/n_products;

    v2 -= v1;

    pcout << "Time for " << n_products <<  " matrix-vector products: MF = "
          << std::setprecision(3) << std::scientific << time_mf
          << "s, SpMV = "
          << std::setprecision(3) << std::scientific << time_sparse
          << "s, verify = "
          << std::setprecision(2) << std::scientific
          << v2.linfty_norm()/v1.linfty_norm()
          << std::endl;
    std::cout << std::setprecision(4);
    std::cout.unsetf(std::ios_base::floatfield);
  }



  template <int dim>
  void LaplaceProblem<dim>::run ()
  {
    for (unsigned int cycle=0; cycle<9-dim; ++cycle)
      {
        {
          TimerOutput::Scope t(timer_output, "Create/refine mesh");
          if (cycle == 0)
            {
              GridGenerator::hyper_cube (triangulation, 0., 1.);
              triangulation.refine_global (3-dim);
            }
          triangulation.refine_global (1);
        }

        setup_system ();
        assemble_rhs ();
        solve ();
        //compare_performance_matrix_vector();
      };
  }
}




int main (int argc, char *argv[])
{
  try
    {
      using namespace Step37;

      Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

      LaplaceProblem<dimension> laplace_problem;
      laplace_problem.run ();
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
