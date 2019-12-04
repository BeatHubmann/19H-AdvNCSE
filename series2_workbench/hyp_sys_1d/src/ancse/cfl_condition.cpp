#include <ancse/cfl_condition.hpp>

#include <Eigen/Dense>


// define FVM CFL condition functions defined in cfl_condition.hpp
FVMCFLCondition::FVMCFLCondition(const Grid& grid,
                                 const std::shared_ptr<Model>& model,
                                 double cfl_number)
    : grid(grid), model(model), cfl_number(cfl_number) {}

double FVMCFLCondition::operator() (const Eigen::MatrixXd& u) const
{
    const int n_cells= grid.n_cells;
    const int n_ghost= grid.n_ghost;

    const double dx= grid.dx;

    double a_max= 0.0;

    for (int i= n_ghost; i < n_cells - n_ghost; ++i)
        a_max= std::max(a_max, model->max_eigenvalue(u.col(i)));
    
    return cfl_number * dx / a_max;
}

// define DG CFL condition functions defined in cfl_condition.hpp
DGCFLCondition::DGCFLCondition(const Grid& grid,
                               const std::shared_ptr<Model>& model,
                               const DGHandler& dg_handler,
                               double cfl_number)
    : grid(grid), model(model), dg_handler(dg_handler), cfl_number(cfl_number) {}

double DGCFLCondition::operator() (const Eigen::MatrixXd& u) const
{
    const int n_cells= grid.n_cells;
    const int n_ghost= grid.n_ghost;

    const double dx= grid.dx;

    double a_max= 0.0;

    Eigen::MatrixXd u0= dg_handler.build_cell_avg(u);

    for (int i= n_ghost; i < n_cells - n_ghost; ++i)
        a_max= std::max(a_max, model->max_eigenvalue(u0.col(i)));
    
    return cfl_number * dx / a_max;
}

/// make CFL condition for FVM
std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid& grid,
                   const std::shared_ptr<Model>& model,
                   double cfl_number)
{
    // implement this 'factory' for your CFL condition.
    return std::make_shared<FVMCFLCondition>(grid, model, cfl_number);
}

/// make CFL condition for DG
std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid& grid,
                   const std::shared_ptr<Model>& model,
                   const DGHandler& dg_handler,
                   double cfl_number)
{
    // implement this 'factory' for your CFL condition.
    return std::make_shared<DGCFLCondition>(grid, model, dg_handler, cfl_number);
}
