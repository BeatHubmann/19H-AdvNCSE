#ifndef FVMSCALAR1D_RUNGE_KUTTA_HPP
#define FVMSCALAR1D_RUNGE_KUTTA_HPP

#include <Eigen/Dense>
#include <map>
#include <memory>

#include <ancse/boundary_condition.hpp>
#include <ancse/rate_of_change.hpp>

/// Interface advancing the solution of a PDE by one step.
class TimeIntegrator {
  public:
    virtual ~TimeIntegrator() = default;

    /// Update 'u1' starting from 'u0' using a time-step `dt`.
    virtual void operator()(Eigen::VectorXd &u1,
                            const Eigen::VectorXd &u0,
                            double dt) const = 0;
};

/// Runge-Kutta time integration methods.
/** Note: not all methods are Runge-Kutta methods. Therefore, 'RungeKutta' can
 *        not be the interface. It should be the base for all (explicit)
 *        Runge-Kutta methods.
 */
class RungeKutta : public TimeIntegrator {
  public:
    RungeKutta(std::shared_ptr<RateOfChange> rate_of_change,
               std::shared_ptr<BoundaryCondition> boundary_condition)
        : rate_of_change(std::move(rate_of_change)),
          boundary_condition(std::move(boundary_condition)) {}

  protected:
    std::shared_ptr<RateOfChange> rate_of_change;
    std::shared_ptr<BoundaryCondition> boundary_condition;
};

/// First order explicit time integration.
class ForwardEuler : public RungeKutta {
  private:
    using super = RungeKutta;

  public:
    ForwardEuler(std::shared_ptr<RateOfChange> rate_of_change,
                 std::shared_ptr<BoundaryCondition> boundary_condition,
                 int n_cells);

    virtual void operator()(Eigen::VectorXd &u1,
                            const Eigen::VectorXd &u0,
                            double dt) const override;

  private:
    mutable Eigen::VectorXd dudt;
};

//// ANCSE_CUT_START_TEMPLATE
//----------------SSP2DeclBegin----------------
class SSP2 : public RungeKutta {
  private:
    using super = RungeKutta;

  public:
    SSP2(std::shared_ptr<RateOfChange> rate_of_change,
         std::shared_ptr<BoundaryCondition> boundary_condition,
         int n_cells);

    virtual void operator()(Eigen::VectorXd &u1,
                            const Eigen::VectorXd &u0,
                            double dt) const override;

  private:
    mutable Eigen::VectorXd u_star;
    mutable Eigen::VectorXd dudt;
};
//----------------SSP2DeclEnd----------------
//// ANCSE_END_TEMPLATE

std::shared_ptr<RungeKutta>
make_runge_kutta(const std::shared_ptr<RateOfChange> &rate_of_change,
                 const std::shared_ptr<BoundaryCondition> &boundary_condition,
                 int n_cells);

#endif // FVMSCALAR1D_RUNGE_KUTTA_HPP
