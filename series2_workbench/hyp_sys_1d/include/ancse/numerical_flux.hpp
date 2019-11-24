#ifndef HYPSYS1D_NUMERICAL_FLUX_HPP
#define HYPSYS1D_NUMERICAL_FLUX_HPP

#include <memory>
#include <iostream>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/simulation_time.hpp>

/// Central flux.
/** This flux works does not depend on the model.
 * It is also unconditionally a bad choice.
 */

// ANCSE lecture notes 4.2.2 (p43): Central schemes:
class CentralFlux {
  public:
    // Note: the interface for creating fluxes will give you access
    //       to the following:
    //         - model
    //         - grid
    //         - shared_ptr to simulation_time
    //       Therefore, try to only use a subset of those three in your
    //       constructors.
    explicit CentralFlux(const std::shared_ptr<Model> &model)
        : model(model) {}

    /// Compute the numerical flux given the left and right trace.
    Eigen::VectorXd operator()(const Eigen::VectorXd &uL,
                               const Eigen::VectorXd &uR) const
    {
        auto fL = model->flux(uL);
        auto fR = model->flux(uR);

        return 0.5 * (fL + fR);
    }

  private:
    std::shared_ptr<Model> model;
};

// ANCSE lecture notes 4.2.4 (p44): Rusanov scheme (1961):
class Rusanov
{
    public:
        explicit Rusanov(const std::shared_ptr<Model>& model) : model(model) {}

        Eigen::VectorXd operator()(const Eigen::VectorXd& uL, const Eigen::VectorXd& uR) const
        {
            const double a= std::max(model->max_eigenvalue(uL), model->max_eigenvalue(uR));
            auto fL= model->flux(uL);
            auto fR= model->flux(uR);

            return 0.5 * ((fL + fR) - a * (uR - uL));
        }

    private:
        std::shared_ptr<Model> model;
};

// ANCSE lecture notes 4.2.3 (p44): Lax-Friedrichs scheme:
class LxF
{
    public:
        explicit LxF(const std::shared_ptr<Model>& model,
                     const Grid& grid,
                     std::shared_ptr<SimulationTime> simulation_time)
                   : model(model),
                     grid(grid),
                     simulation_time(std::move(simulation_time)) {}

        Eigen::VectorXd operator()(const Eigen::VectorXd& uL, const Eigen::VectorXd& uR) const
        {
            const double dx= grid.dx;
            const double dt= simulation_time->dt;
            const double a= dx / dt;
            auto fL= model->flux(uL);
            auto fR= model->flux(uR);

            return 0.5 * ((fL + fR) - a * (uR - uL));
        }

    private:
        std::shared_ptr<Model> model;
        Grid grid;
        std::shared_ptr<SimulationTime> simulation_time;
};

// ANCSE lecture notes 9.2.1 (p118): Roe's scheme (1981):
class Roe
{
    public:
        explicit Roe(const std::shared_ptr<Model>& model) : model(model) {}

        Eigen::VectorXd operator()(const Eigen::VectorXd uL, const Eigen::VectorXd uR) const
        {
        
            auto fL= model->flux(uL);
            auto fR= model->flux(uR);
            // TODO: Calculate Roe matrix/flux
            return fL;
        }
    
    private:
        std::shared_ptr<Model> model;
};

class HLL
{

};

class HLLC
{

};

#endif // HYPSYS1D_NUMERICAL_FLUX_HPP
