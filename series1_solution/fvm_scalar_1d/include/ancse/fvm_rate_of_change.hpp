#ifndef FVMSCALAR1D_FVM_RATE_OF_CHANGE_HPP
#define FVMSCALAR1D_FVM_RATE_OF_CHANGE_HPP

#include <Eigen/Dense>
#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/rate_of_change.hpp>
#include <ancse/simulation_time.hpp>
#include <memory>

//----------------FVMRateOfChangeBegin----------------
/// Compute the rate of change due to FVM.
/** The semidiscrete approximation of a PDE using FVM is
 *      du_i/dt = - (F_{i+0.5} - F_{i-0.5}) / dx.
 *  This computes the right hand side of the ODE.
 *
 * @tparam NumericalFlux see e.g. `CentralFlux`.
 * @tparam Reconstruction see e.g. `PWConstantReconstruction`.
 */
template <class NumericalFlux, class Reconstruction>
class FVMRateOfChange : public RateOfChange {
  public:
    FVMRateOfChange(const Grid &grid,
                    const NumericalFlux &numerical_flux,
                    const Reconstruction &reconstruction)
        : grid(grid),
          numerical_flux(numerical_flux),
          reconstruction(reconstruction) {}

    virtual void operator()(Eigen::VectorXd &dudt,
                            const Eigen::VectorXd &u0) const override {
        //// ANCSE_CUT_START_TEMPLATE
        //// ANCSE_COMMENT implement the flux loop here.
        auto n_cells = grid.n_cells;
        auto n_ghost = grid.n_ghost;

        double dx = grid.dx;
        double fL, fR = 0.0;
        double uL, uR;

        for (int i = n_ghost - 1; i < n_cells - n_ghost; ++i) {
            std::tie(uL, uR) = reconstruction(u0, i);

            fL = fR;
            fR = numerical_flux(uL, uR);

            dudt[i] = (fL - fR) / dx;
        }
        //// ANCSE_END_TEMPLATE
    }

  private:
    Grid grid;
    NumericalFlux numerical_flux;
    Reconstruction reconstruction;
};
//----------------FVMRateOfChangeEnd----------------

std::shared_ptr<RateOfChange>
make_fvm_rate_of_change(const Grid &grid,
                        const Model &model,
                        const std::shared_ptr<SimulationTime> &simulation_time);

#endif // FVMSCALAR1D_FVM_RATE_OF_CHANGE_HPP
