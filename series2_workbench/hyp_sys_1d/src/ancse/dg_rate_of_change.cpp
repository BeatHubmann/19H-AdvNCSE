#include <ancse/dg_rate_of_change.hpp>

#include <Eigen/Dense>
#include <ancse/config.hpp>
#include <ancse/polynomial_basis.hpp>
#include <ancse/dg_handler.hpp>
#include <ancse/numerical_flux.hpp>
#include <fmt/format.h>


/// DG numerical flux term
template <class NumericalFlux>
void DGRateOfChange<NumericalFlux>
:: eval_numerical_flux (Eigen::MatrixXd &dudt,
                        const Eigen::MatrixXd &u0) const
{
    // implement the loop for DG numerical flux term.
    const int n_cells= grid.n_cells;
    const int n_ghost= grid.n_ghost;

    const int n_vars= model->get_nvars();

    const double dx= grid.dx;
    Eigen::VectorXd fL= Eigen::VectorXd::Zero(n_vars), fR= Eigen::VectorXd::Zero(n_vars);
    Eigen::VectorXd uL, uR;

    for (int i= n_ghost - 1; i < n_cells - n_ghost; ++i)
    {
        uL= u0.col(i);
        uR= u0.col(i + 1);
        
        fL= fR;
        fR= numerical_flux(uL, uR);

        dudt.col(i)= (fL - fR) / dx;
    }
}

/// DG volume integral term
template <class NumericalFlux>
void DGRateOfChange<NumericalFlux>
:: eval_volume_integral(Eigen::MatrixXd &dudt,
                        const Eigen::MatrixXd &u0) const
{
    // implement the loop for DG volume integral.
}


#define REGISTER_NUMERICAL_FLUX(token, FluxType, flux)          \
    if (config["flux"] == (token)) {                            \
        return std::make_shared< DGRateOfChange<FluxType> >(    \
            grid, model, flux, poly_basis, dg_handler);                     \
    }


std::shared_ptr<RateOfChange> make_dg_rate_of_change(
    const nlohmann::json &config,
    const Grid &grid,
    const std::shared_ptr<Model> &model,
    const PolynomialBasis &poly_basis,
    const DGHandler &dg_handler,
    const std::shared_ptr<SimulationTime> &simulation_time)
{
    // Register the other numerical fluxes.

    REGISTER_NUMERICAL_FLUX("central_flux", CentralFlux, CentralFlux(model))
    REGISTER_NUMERICAL_FLUX("rusanov", Rusanov, Rusanov(model))
    REGISTER_NUMERICAL_FLUX("lax_friedrichs",
                            LaxFriedrichs,
                            LaxFriedrichs(grid, model, simulation_time))
    REGISTER_NUMERICAL_FLUX("roe", Roe, Roe(model))
    REGISTER_NUMERICAL_FLUX("hll", HLL, HLL(model))
    REGISTER_NUMERICAL_FLUX("hllc", HLLc, HLLc(model))

    throw std::runtime_error(
        fmt::format("Unknown numerical flux. {}",
                    std::string(config["flux"])));
}

#undef REGISTER_NUMERICAL_FLUX
