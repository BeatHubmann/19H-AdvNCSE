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


/// Rusanov flux.
/**
 * 
 */
class RusanovFlux
{
    public:
        explicit RusanovFlux(const Model &model) : model(model) {}

        double operator()(double uL, double uR) const
        {
            auto fL= model.flux(uL);
            auto fR= model.flux(uR);
            auto s= std::max(std::abs(model.max_eigenvalue(uL)),
                                 std::abs(model.max_eigenvalue(uR)));            
            return 0.5 * ((fL + fR) - s * (uR - uL));
        }

    private:
        Model model;
};

/// Lax Friedrichs flux.
/**
 * 
 */
class LxFFlux
{
    public:
        explicit LxFFlux(const Model &model, const double dx, const double dt) : model(model), s(dx/dt) {}

        double operator()(double uL, double uR) const
        {
            auto fL= model.flux(uL);
            auto fR= model.flux(uR);
            return 0.5 * ((fL + fR) - s * (uR - uL));
        }

    private:
        Model model;
        double s;
};

#endif // FVMSCALAR1D_NUMERICAL_FLUX_HPP
