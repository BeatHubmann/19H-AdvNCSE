#ifndef FVMSCALAR1D_NUMERICAL_FLUX_HPP
#define FVMSCALAR1D_NUMERICAL_FLUX_HPP

#include <memory>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/simulation_time.hpp>

/// Central flux.
/** This flux works does not depend on the model. It is also unconditionally a
 * bad choice.
 */
class CentralFlux {
  public:
    // Note: the interface for creating fluxes will give you access to the
    //       following:
    //         - model
    //         - grid
    //         - shared_ptr to simulation_time
    //       Therefore, try to only use a subset of those three in your
    //       constructors.
    explicit CentralFlux(const Model &model) : model(model) {}

    /// Compute the numerical flux given the left and right trace.
    double operator()(double uL, double uR) const {
        auto fL = model.flux(uL);
        auto fR = model.flux(uR);

        return 0.5 * (fL + fR);
    }

  private:
    Model model;
};

//// ANCSE_CUT_START_TEMPLATE
//----------------FluxLFBegin----------------
/// Lax-Friedrichs numerical flux.
/** This flux works for any model. */
class LaxFriedrichs {
  public:
    // Note: This version is a bit tricky. A numerical flux should be
    //       a function of the two trace values at the interface, i.e. what we
    //       call `uL`, `uR`. However, it requires 'dt' and 'dx'. Therefore,
    //       these need to be made available to the flux. This is one of the
    //       reasons why `SimulationTime`.
    LaxFriedrichs(const Grid &grid,
                  const Model &model,
                  std::shared_ptr<SimulationTime> simulation_time)
        : simulation_time(std::move(simulation_time)),
          grid(grid),
          model(model) {}

    double operator()(double uL, double uR) const {
        double dx = grid.dx;
        double dt = simulation_time->dt;

        auto fL = model.flux(uL);
        auto fR = model.flux(uR);

        return 0.5 * ((fL + fR) - dx / dt * (uR - uL));
    }

  private:
    std::shared_ptr<SimulationTime> simulation_time;
    Grid grid;
    Model model;
};
//----------------FluxLFEnd----------------
//// ANCSE_CUT_END_TEMPLATE

//// ANCSE_CUT_START_TEMPLATE
//----------------FluxRusanovBegin----------------
/// Rusanov's flux (or local Lax-Friedrichs).
/** This flux works for any model. */
class Rusanov {
  public:
    explicit Rusanov(const Model &model) : model(model) {}

    double operator()(double uL, double uR) const {
        double a = std::max(model.max_eigenvalue(uL), model.max_eigenvalue(uR));

        auto fL = model.flux(uL);
        auto fR = model.flux(uR);

        return 0.5 * ((fL + fR) - a * (uR - uL));
    }

  private:
    Model model;
};
//----------------FluxRusanovEnd----------------
//// ANCSE_CUT_END_TEMPLATE

//// ANCSE_CUT_START_TEMPLATE
//----------------FluxRoeBegin----------------
/// The infamous Roe flux.
/** Works for any scalar conservation law.
 */
class Roe {
  public:
    explicit Roe(const Model &model) : model(model) {}

    double operator()(double uL, double uR) const {
        auto fL = model.flux(uL);
        auto fR = model.flux(uR);

        bool is_equal = std::abs(uL - uR) < atol;
        double a = is_equal ? model.max_eigenvalue(uL) : (fL - fR) / (uL - uR);
        // The value of a doesn't matter in the case is_equal

        return (a >= 0.0 ? fL : fR);
    }

  private:
    Model model;
    double atol = 1e-10;
};
//----------------FluxRoeBegin----------------
//// ANCSE_END_TEMPLATE

//// ANCSE_CUT_START_TEMPLATE
//----------------FluxEOBegin----------------
/// The numerical flux by Engquist and Osher.
/** Engquist-Osher requires the some pen&paper computations specific to the
 *  convervation law at hand.
 *
 *  This version is for the Burgers' equation.
 */
class EngquistOsher {
  public:
    explicit EngquistOsher(const Model &model) : model(model) {}

    double operator()(double uL, double uR) const {
        double omega = 0.0;
        double f_plus = model.flux(std::max(uL, omega));
        double f_minus = model.flux(std::min(uR, omega));

        return f_plus + f_minus;
    }

  private:
    Model model;
};
//----------------FluxEOEnd----------------
//// ANCSE_END_TEMPLATE

//// ANCSE_CUT_START_TEMPLATE
//----------------FluxGodunovBegin----------------
/// The numerical flux due to Godunov
/** Godunov's flux requires the some pen&paper computations specific to the
 *  convervation law at hand.
 *
 *  This version is for the Burgers' equation.
 */
class Godunov {
  public:
    explicit Godunov(const Model &model) : model(model) {}

    double operator()(double uL, double uR) const {
        double omega = 0.0;

        double f_plus = model.flux(std::max(uL, omega));
        double f_minus = model.flux(std::min(uR, omega));

        return std::max(f_plus, f_minus);
    }

  private:
    Model model;
};
//----------------FluxGodunovEnd----------------
//// ANCSE_END_TEMPLATE

#endif // FVMSCALAR1D_NUMERICAL_FLUX_HPP
