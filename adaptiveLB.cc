#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
 
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/affine_constraints.h> //from step-6
 
#include <deal.II/grid/tria.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h> //from step-6
#include <deal.II/grid/grid_refinement.h> //from step-6

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
 
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>
 
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h> //from step-6
 
#include <fstream>
#include <iostream>

namespace adaptiveLB
{
  using namespace dealii;
 
 

  template <int spacedim>
  class adaptiveLBProblem
  {
  public:
    adaptiveLBProblem(const unsigned degree = 2);
    void run();
 
  private:
    static constexpr unsigned int dim = spacedim - 1;
 
    void make_grid_and_dofs(); // in step-6 is called setup_system
    void assemble_system();
    void solve();
    void output_results(const unsigned int cycle) const; // we embed the cycle number in the output file name as in step-6
    void compute_error() const;
    void refine_grid(); // add this function as in step-6
    
 
    Triangulation<dim, spacedim> triangulation;
    FE_Q<dim, spacedim>          fe;
    DoFHandler<dim, spacedim>    dof_handler;
    MappingQ<dim, spacedim>      mapping;

    AffineConstraints<double> constraints; // constraint for hanging nodes
 
    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
 
    Vector<double> solution;
    Vector<double> system_rhs;
  };
 
// maybe this later, go to make grid and dofs
 
  template <int dim>
  class Solution : public Function<dim>
  {
  public:
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
 
    virtual Tensor<1, dim>
    gradient(const Point<dim> & p,
             const unsigned int component = 0) const override;
  };
 
 
  template <>
  double Solution<2>::value(const Point<2> &p, const unsigned int) const
  {
    return (-2. * p(0) * p(1));
  }
 
 
  template <>
  Tensor<1, 2> Solution<2>::gradient(const Point<2> &p,
                                     const unsigned int) const
  {
    Tensor<1, 2> return_value;
    return_value[0] = -2. * p(1) * (1 - 2. * p(0) * p(0));
    return_value[1] = -2. * p(0) * (1 - 2. * p(1) * p(1));
 
    return return_value;
  }
 
 
  template <>
  double Solution<3>::value(const Point<3> &p, const unsigned int) const
  {
    return (std::sin(numbers::PI * p(0)) * std::cos(numbers::PI * p(1)) *
            exp(p(2)));
  }
 
 
  template <>
  Tensor<1, 3> Solution<3>::gradient(const Point<3> &p,
                                     const unsigned int) const
  {
    using numbers::PI;
 
    Tensor<1, 3> return_value;
 
    return_value[0] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
    return_value[1] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
    return_value[2] = sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 
    return return_value;
  }
 
 
 
  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    virtual double value(const Point<dim> & p,
                         const unsigned int component = 0) const override;
  };
 
  template <>
  double RightHandSide<2>::value(const Point<2> &p,
                                 const unsigned int /*component*/) const
  {
    return (-8. * p(0) * p(1));
  }
 
 
  template <>
  double RightHandSide<3>::value(const Point<3> &p,
                                 const unsigned int /*component*/) const
  {
    using numbers::PI;
 
    Tensor<2, 3> hessian;
 
    hessian[0][0] = -PI * PI * sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
    hessian[1][1] = -PI * PI * sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
    hessian[2][2] = sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 
    hessian[0][1] = -PI * PI * cos(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
    hessian[1][0] = -PI * PI * cos(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
 
    hessian[0][2] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
    hessian[2][0] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 
    hessian[1][2] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
    hessian[2][1] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
 
    Tensor<1, 3> gradient;
    gradient[0] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
    gradient[1] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2));
    gradient[2] = sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2));
 
    Point<3> normal = p;
    normal /= p.norm();
 
    return (-trace(hessian) + 2 * (gradient * normal) +
            (hessian * normal) * normal);
  }
 
 
 
  template <int spacedim>
  adaptiveLBProblem<spacedim>::adaptiveLBProblem(
    const unsigned degree)
    : fe(degree)
    , dof_handler(triangulation)
    , mapping(degree)
  {}
 
 
  template <int spacedim>
  void adaptiveLBProblem<spacedim>::make_grid_and_dofs()
  {
    dof_handler.distribute_dofs(fe);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.clear();

    DoFTools::make_hanging_node_constraints(dof_handler, constraints);

    // Note that we do not to apply the boundary conditions after assembly,
    // instead we put all constraints on our function space in the AffineConstraints object.
    VectorTools::interpolate_boundary_values(dof_handler,
                                           0,
                                           Functions::ZeroFunction<spacedim>(), // changed dim in spacedim
                                           constraints);

    constraints.close();

    // I think that this part will go in the run function, as in the step-6
    /*
    {
      Triangulation<spacedim> volume_mesh;
      GridGenerator::half_hyper_ball(volume_mesh);
 
      const std::set<types::boundary_id> boundary_ids = {0};
 
      GridGenerator::extract_boundary_mesh(volume_mesh,
                                           triangulation,
                                           boundary_ids);
    }
    triangulation.set_all_manifold_ids(0);
    triangulation.set_manifold(0, SphericalManifold<dim, spacedim>());
 
    triangulation.refine_global(4);
    */
 
    std::cout << "Surface mesh has " << triangulation.n_active_cells()
              << " cells." << std::endl;
 
 
    std::cout << "Surface mesh has " << dof_handler.n_dofs()
              << " degrees of freedom." << std::endl;
 
    DynamicSparsityPattern dsp(dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler,
                                  dsp,
                                  constraints,
                                  /*keep_constrained_dofs = */ false);
    sparsity_pattern.copy_from(dsp);
 
    system_matrix.reinit(sparsity_pattern);
  }
 
 
 
  template <int spacedim>
  void adaptiveLBProblem<spacedim>::assemble_system()
  {
    system_matrix = 0;
    system_rhs    = 0;
 
    const QGauss<dim>       quadrature_formula(2 * fe.degree);
    FEValues<dim, spacedim> fe_values(mapping,
                                      fe,
                                      quadrature_formula,
                                      update_values | update_gradients |
                                        update_quadrature_points |
                                        update_JxW_values);
 
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
    const unsigned int n_q_points    = quadrature_formula.size();
 
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);
 
    std::vector<double>                  rhs_values(n_q_points);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 
    RightHandSide<spacedim> rhs;
 
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        cell_matrix = 0;
        cell_rhs    = 0;
 
        fe_values.reinit(cell);
 
        rhs.value_list(fe_values.get_quadrature_points(), rhs_values);
 
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = 0; j < dofs_per_cell; ++j)
            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
              cell_matrix(i, j) += fe_values.shape_grad(i, q_point) *
                                   fe_values.shape_grad(j, q_point) *
                                   fe_values.JxW(q_point);
 
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            cell_rhs(i) += fe_values.shape_value(i, q_point) *
                           rhs_values[q_point] * fe_values.JxW(q_point);
 
        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);
        // not needed more, we have do it in the line above
        /* 
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              system_matrix.add(local_dof_indices[i],
                                local_dof_indices[j],
                                cell_matrix(i, j));
 
            system_rhs(local_dof_indices[i]) += cell_rhs(i);
          }
        */
      }
    //this part is not needed anymore, we have put it in the constraints
    /*
    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(
      mapping, dof_handler, 0, Solution<spacedim>(), boundary_values);
 
    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix, solution, system_rhs, false);
    */
  }
 
 
 
 
  template <int spacedim>
  void adaptiveLBProblem<spacedim>::solve()
  {
    SolverControl solver_control(solution.size(), 1e-7 * system_rhs.l2_norm());
    SolverCG<Vector<double>> cg(solver_control);
 
    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);
 
    cg.solve(system_matrix, solution, system_rhs, preconditioner);
    
    constraints.distribute(solution);
  }
 
  template <int spacedim> // changed dim in spacedim
  void adaptiveLBProblem<spacedim>::refine_grid() // changed dim in spacedim
  {
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
  
    KellyErrorEstimator<dim,spacedim>::estimate(dof_handler, // added space dim in the KellyErrorEstimator, dont know if is correct. Maybe remove it, look at step -6.
                                      QGauss<dim - 1>(fe.degree + 1), // dim is now spacedim-1, so its correct dim-1
                                      {},
                                      solution,
                                      estimated_error_per_cell);
  
    GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                    estimated_error_per_cell,
                                                    0.3,
                                                    0.03);
  
    triangulation.execute_coarsening_and_refinement();
  }
 
 
  template <int spacedim>
  void adaptiveLBProblem<spacedim>::output_results(const unsigned int cycle) const
  {
    {
      // its a little bit tricky, to me, to understand if this part is correct or not. Maybe read better the step-6.
      GridOut               grid_out;
      std::ofstream         output("grid-" + std::to_string(cycle) + ".gnuplot");
      GridOutFlags::Gnuplot gnuplot_flags(false, 5);
      grid_out.set_flags(gnuplot_flags);
      grid_out.write_gnuplot(triangulation, output, &mapping);
    }

    {
      DataOut<dim, spacedim> data_out;
      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution,
                              "solution",
                              DataOut<dim, spacedim>::type_dof_data);
      data_out.build_patches(mapping, mapping.get_degree());
  
      const std::string filename =
        "solution" + std::to_string(spacedim) + "d-cycle_" + std::to_string(cycle)+ ".vtk";
      std::ofstream output(filename);
      data_out.write_vtk(output);
    }
  }
 
 
 
 
  template <int spacedim> // dont know if i want to put it in this program.
  void adaptiveLBProblem<spacedim>::compute_error() const
  {
    Vector<float> difference_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(mapping,
                                      dof_handler,
                                      solution,
                                      Solution<spacedim>(),
                                      difference_per_cell,
                                      QGauss<dim>(2 * fe.degree + 1),
                                      VectorTools::H1_norm);
 
    double h1_error = VectorTools::compute_global_error(triangulation,
                                                        difference_per_cell,
                                                        VectorTools::H1_norm);
    std::cout << "H1 error = " << h1_error << std::endl;
  }
 
 
 
 
  template <int spacedim>
  void adaptiveLBProblem<spacedim>::run()
  {
    for (unsigned int cycle = 0; cycle < 8; ++cycle)
    {
      std::cout << "Cycle " << cycle << ':' << std::endl;
 
      if (cycle == 0)
        {
          {
            Triangulation<spacedim> volume_mesh;
            GridGenerator::half_hyper_ball(volume_mesh);
      
            const std::set<types::boundary_id> boundary_ids = {0};
      
            GridGenerator::extract_boundary_mesh(volume_mesh,
                                                triangulation,
                                                boundary_ids);
          }
          triangulation.set_all_manifold_ids(0);
          triangulation.set_manifold(0, SphericalManifold<dim, spacedim>());
      
          triangulation.refine_global(1);
          // the following is not needed anymore
          /*
          GridGenerator::hyper_ball(triangulation);
          triangulation.refine_global(1);
          */
        }
      else
        refine_grid();
 
 
      std::cout << "   Number of active cells:       "
                << triangulation.n_active_cells() << std::endl;
 
      make_grid_and_dofs();
 
      std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs()
                << std::endl;
 
      assemble_system();
      solve();
      output_results(cycle);
      compute_error(); // dont know if i want to put it in this program.
    }

  }
} // namespace adaptiveLB
 
 
 
int main()
{
  try
    {
      using namespace adaptiveLB;
 
      adaptiveLBProblem<3> laplace_beltrami;
      laplace_beltrami.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
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
      std::cerr << std::endl
                << std::endl
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