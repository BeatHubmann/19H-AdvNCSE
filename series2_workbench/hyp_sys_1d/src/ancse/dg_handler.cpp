#include <ancse/dg_handler.hpp>


/// build solution from DG coefficients and the basis
/// pre-evaluated at a certain point
Eigen::VectorXd DGHandler
:: build_sol(const Eigen::VectorXd& u,
             const Eigen::VectorXd& basis) const
{
    Eigen::VectorXd uSol(n_vars);

    // looping over m(=3 for 1D Euler) n_vars in the coefficient column vector with stride size n_coeff:
    for (int l= 0; l < n_vars; ++l)
        uSol(l)= u.segment(l * n_coeff, n_coeff).dot(basis);

    return uSol;
}

/// build solution from DG coefficients at a given reference point
Eigen::VectorXd DGHandler 
:: build_sol(const Eigen::VectorXd& u,
             double x) const
{
    // convert physical point x in cell K_i into reference point xi:
    const double xi= reference_point(grid, x);
    
    return build_sol(u, poly_basis(xi));
}

/// build cell average
Eigen::MatrixXd DGHandler
:: build_cell_avg (const Eigen::MatrixXd& u) const
{
    auto n_cells = u.cols();
    Eigen::MatrixXd u0= Eigen::MatrixXd::Zero(n_vars, n_cells);

    // iterating over matrix u0 and filling from matrix u:
    for (int i= 0; i < n_vars; ++i)
        for (int j= 0; j < n_cells; ++j)
            u0(i, j)+= u(i * n_coeff, j);

    return u0 / grid.dx;
}

/// build split solution uSol_m = u0 + um, uSol_p = u0 - up
/// from DG coefficients
std::tuple <Eigen::MatrixXd, Eigen::MatrixXd, Eigen::MatrixXd>
DGHandler :: build_split_sol(const Eigen::MatrixXd& u) const
{
    auto n_cells = u.cols();
    if (n_coeff > 3) {
        throw std::runtime_error(
            "Limiter not implemented for higher than 3rd order");
    }
    
    auto u0 = build_cell_avg(u);
    Eigen::MatrixXd um = Eigen::MatrixXd::Zero (n_vars, n_cells);
    Eigen::MatrixXd up = Eigen::MatrixXd::Zero (n_vars, n_cells);


    return {std::move(u0), std::move(um), std::move(up)};
}

/// build DG coefficients from uSol_m = u0 + um, uSol_p = u0 - up
void DGHandler :: compute_limit_coeffs (Eigen::MatrixXd &u,
                                        Eigen::MatrixXd &um,
                                        Eigen::MatrixXd &up) const
{
    if (n_coeff == 1) {
        return;
    }
    else if (n_coeff > 3) {
        throw std::runtime_error(
            "Limiter not implemented for higher than 3rd order");
    }
    
}
