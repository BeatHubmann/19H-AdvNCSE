#ifndef HYPSYS1D_FVM_RATE_OF_CHANGE_HPP
#define HYPSYS1D_FVM_RATE_OF_CHANGE_HPP

#include <memory>
#include <Eigen/Dense>

#include <ancse/config.hpp>
#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/rate_of_change.hpp>
#include <ancse/simulation_time.hpp>

/// Compute the rate of change due to FVM.
/** The semidiscrete approximation of a PDE using FVM is
 *      du_i/dt = - (F_{i+0.5} - F_{i-0.5}) / dx.
 *  This computes the right hand side of the ODE.
 *
 * @tparam NumericalFlux see e.g. `CentralFlux`.
 * @tparam Reconstruction see e.g. `PWConstantReconstruction`.
 */
template <class NumericalFlux, class Reconstruction>
class FVMRateOfChange : public RateOfChange
{
    public:
        FVMRateOfChange(const Grid& grid,
                        const std::shared_ptr<Model>& model,
                        const NumericalFlux& numerical_flux,
                        const Reconstruction& reconstruction)
            : grid(grid),
              model(model),
              numerical_flux(numerical_flux),
              reconstruction(reconstruction) {}

        virtual void operator()(Eigen::MatrixXd& dudt,
                                const Eigen::MatrixXd& u0) const override
        {
            // implement the flux loop here.
            const int n_cells= grid.n_cells;
            const int n_ghost= grid.n_ghost;

            const int n_vars= model->get_nvars();

            const double dx= grid.dx;
            Eigen::VectorXd fL= Eigen::VectorXd::Zero(n_vars), fR= Eigen::VectorXd::Zero(n_vars);
            Eigen::VectorXd uL, uR;

            for (int i= n_ghost - 1; i < n_cells - n_ghost; ++i)
            {
                std::tie(uL, uR)= reconstruction(u0, i);
             
                fL= fR;
                fR= numerical_flux(uL, uR);

                dudt.col(i)= (fL - fR) / dx;
            }
        }

    private:
        Grid grid;
        std::shared_ptr<Model> model;
        NumericalFlux numerical_flux;
        Reconstruction reconstruction;
};

std::shared_ptr<RateOfChange>
make_fvm_rate_of_change(const nlohmann::json &config,
                        const Grid &grid,
                        const std::shared_ptr<Model> &model,
                        const std::shared_ptr<SimulationTime> &simulation_time);

#endif // HYPSYS1D_FVM_RATE_OF_CHANGE_HPP
