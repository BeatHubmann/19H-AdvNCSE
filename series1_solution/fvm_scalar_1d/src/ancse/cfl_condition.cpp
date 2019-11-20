#include <ancse/cfl_condition.hpp>

//// ANCSE_CUT_START_TEMPLATE
//----------------StandardCFLConditionDefnBegin----------------
StandardCFLCondition::StandardCFLCondition(const Grid &grid,
                                           const Model &model,
                                           double cfl_number)
    : grid(grid), model(model), cfl_number(cfl_number) {}

double StandardCFLCondition::operator()(const Eigen::VectorXd &u) const {

    auto n_cells = grid.n_cells;
    auto n_ghost = grid.n_ghost;

    double a_max = 0.0;
    for (int i = grid.n_ghost; i < n_cells - n_ghost; ++i) {
        a_max = std::max(a_max, model.max_eigenvalue(u[i]));
    }

    return cfl_number * grid.dx / a_max;
}
//----------------StandardCFLConditionDefnEnd----------------
//// ANCSE_END_TEMPLATE

std::shared_ptr<CFLCondition>
make_cfl_condition(const Grid &grid, const Model &model, double cfl_number) {
    //// ANCSE_CUT_START_TEMPLATE
    //// ANCSE_COMMENT implement this 'factory' for your CFL condition.
    return std::make_shared<StandardCFLCondition>(grid, model, cfl_number);
    //// ANCSE_RETURN_VALUE nullptr;
    //// ANCSE_END_TEMPLATE
}
