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

        Eigen::VectorXd operator() (const Eigen::VectorXd& uL, const Eigen::VectorXd& uR) const
        {
            const double a= std::max(model->max_eigenvalue(uL), model->max_eigenvalue(uR));
            auto fL= model->flux(uL);
            auto fR= model->flux(uR);

            return 0.5 * (  (fL + fR)
                          - a * (uR - uL));
        }

    private:
        std::shared_ptr<Model> model;
};

// ANCSE lecture notes 4.2.3 (p44): Lax-Friedrichs scheme:
class LaxFriedrichs
{
    public:
        explicit LaxFriedrichs(const Grid& grid,
            const std::shared_ptr<Model>& model,
            std::shared_ptr<SimulationTime> simulation_time)
            : grid(grid),
              model(model),
              simulation_time(std::move(simulation_time)) {}

        Eigen::VectorXd operator() (const Eigen::VectorXd& uL, const Eigen::VectorXd& uR) const
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

// ANCSE lecture notes 9.2.1 (p118): Roe's scheme (1981) cf. Toro p.357:
class Roe
{
    public:
        explicit Roe(const std::shared_ptr<Model>& model) : model(model) {}

        Eigen::VectorXd operator() (const Eigen::VectorXd uL, const Eigen::VectorXd uR) const
        {
        
            auto fL= model->flux(uL);
            auto fR= model->flux(uR);

            auto gamma= model->get_gamma();
            auto rhoL= model->rho(uL);
            auto rhoR= model->rho(uR);
            auto vL= model->v(uL);
            auto vR= model->v(uR);
            auto HL= model->H(uL);
            auto HR= model->H(uR);

            const double rhoL_sqrt= std::sqrt(rhoL);
            const double rhoR_sqrt= std::sqrt(rhoR);

            const int n_vars= model->get_nvars();

            // Step 1: Calculate Roe average values:
            const double rho_avg= 0.5 * (rhoL + rhoR);
            const double v_avg= (rhoL_sqrt * vL + rhoR_sqrt * vR) / (rhoL_sqrt + rhoR_sqrt);
            const double H_avg= (rhoL_sqrt * HL + rhoR_sqrt * HR) / (rhoL_sqrt + rhoR_sqrt);
            const double c_avg_2= (gamma - 1) * (H_avg - 0.5 * v_avg * v_avg);
            const double c_avg= std::sqrt(c_avg_2);

            // Step 2: Calculate averaged eigenvalues:         
            const double lambda_avg_1= v_avg - c_avg;
            const double lambda_avg_2= v_avg;
            const double lambda_avg_5= v_avg + c_avg;

            // Step 3: Compute averaged right eigenvectors:
            Eigen::VectorXd right_eigenvec_1(n_vars), right_eigenvec_2(n_vars), right_eigenvec_5(n_vars); 
            
            right_eigenvec_1 << 1.0,
                                v_avg - c_avg,
                                H_avg - v_avg * c_avg;
            
            right_eigenvec_2 << 1.0,
                                v_avg,
                                0.5 * v_avg * v_avg * v_avg * v_avg; // ARE THESE REALLY v2?
            
            right_eigenvec_5 << 1.0,
                                v_avg + c_avg,
                                H_avg + v_avg * c_avg;



            // Step 4: Compute wave strengths:
            // Step 4a: Need jumps in conserved quantity u_i:
            const double delta_u_1= uR(0) - uL(0);
            const double delta_u_2= uR(1) - uL(1);
            const double delta_u_5= uR(2) - uL(2);
            // Step 4b: Actual wave strengths:
            const double alpha_2= (gamma - 1.0) / c_avg_2 * (delta_u_1 * (H_avg - v_avg * v_avg) + v_avg * delta_u_2 - delta_u_5);
            const double alpha_1= 0.5 / c_avg * (delta_u_1 * (v_avg + c_avg) - delta_u_2 - c_avg * alpha_2);
            const double alpha_5= delta_u_1 - (alpha_1 + alpha_2);

            // Step 5: Assemble above quantities into flux:
            return 0.5 * (  (fL + fR)
                          - alpha_1 * std::abs(lambda_avg_1) * right_eigenvec_1
                          - alpha_2 * std::abs(lambda_avg_2) * right_eigenvec_2
                          - alpha_5 * std::abs(lambda_avg_5) * right_eigenvec_5);
        }
    
    private:
        std::shared_ptr<Model> model;
};

// ANCSE lecture notes 9.3 (p.120): HLL scheme (1983); wave speeds given by (9.27) cf. Toro p.320:
class HLL
{
    public:
        explicit HLL(const std::shared_ptr<Model>& model) : model(model) {}

        Eigen::VectorXd operator() (const Eigen::VectorXd uL, const Eigen::VectorXd uR) const
        {
            auto fL= model->flux(uL);
            auto fR= model->flux(uR);

            // const int n_vars= model->get_nvars();

            // calculate simple choice L and R wave speeds (ANCSE (9.27)):
            double lo_lambdaL, hi_lambdaL, lo_lambdaR, hi_lambdaR;
            std::tie(lo_lambdaL, hi_lambdaL)= model->lo_hi_eigenvalues(uL);
            std::tie(lo_lambdaR, hi_lambdaR)= model->lo_hi_eigenvalues(uR);

            const double sL= std::min(lo_lambdaL, lo_lambdaR);
            const double sR= std::max(hi_lambdaL, hi_lambdaR);

            if (sL >= 0.0)
                return fL;
            else if (sR <= 0.0)
                return fR;
            else
                return (sR * fL - sL * fR + sR * sL * (uR - uL)) / (sR - sL);
        }

    private:
        std::shared_ptr<Model> model;
};

// ANCSE lecture notes 9.4 (p.120): HLLc 3-wave solver (1994) cf. Toro p.322
class HLLc
{
    public:
        explicit HLLc(const std::shared_ptr<Model>& model) : model(model) {}
    
        Eigen::VectorXd operator() (const Eigen::VectorXd uL, const Eigen::VectorXd uR) const
        {
            auto fL= model->flux(uL);
            auto fR= model->flux(uR);

            auto rhoL= model->rho(uL);
            auto rhoR= model->rho(uR);
            auto vL= model->v(uL);
            auto vR= model->v(uR);
            auto pL= model->p(uL);
            auto pR= model->p(uL);

            const int n_vars= model->get_nvars();

            // calculate simple choice L and R wave speeds (ANCSE (9.27)):
            double lo_lambdaL, hi_lambdaL, lo_lambdaR, hi_lambdaR;
            std::tie(lo_lambdaL, hi_lambdaL)= model->lo_hi_eigenvalues(uL);
            std::tie(lo_lambdaR, hi_lambdaR)= model->lo_hi_eigenvalues(uR);

            const double sL= std::min(lo_lambdaL, lo_lambdaR);
            const double sR= std::max(hi_lambdaL, hi_lambdaR);            

            // calculate intermediate estimated M wave speed (ANCSE (9.33), Toro (10.37)):
            const double sM= (pR - pL + rhoL * vL * (sL - vL) - rhoR * vR * (sR - vR)
                             / (rhoL * (sL - vL) - rhoR * (sR - vR)));
            // std::cout << sL << " " << sM << " " << sR << "\n";

            // calculate intermediate pressure p* := p*R = p*R; ANCSE (9.34), Toro (10.30):
            const double pM= pR + rhoR * (sR - vR) * (sM - vR);            

            // setup helper vector D* for calculating intermediate fluxes; Toro (10.40):
            Eigen::VectorXd d(n_vars);
            d << 0.0,
                 1.0,
                 sM;

            // return choice of flux depeding on signal speeds; Toro (10.26), (10.41):
            if (sL >= 0.0)
                return fL;
            else if (sL < 0.0 && sM >= 0.0)
            {
                Eigen::VectorXd uML= (sL * uL - fL + pM * d) / (sL - sM);
                return (fL + sL * (uML - uL));
            }
            else if (sM < 0.0 && sR >= 0.0)
            {
                Eigen::VectorXd uMR= (sR * uR - fR + pM * d) / (sR - sM);
                return (fR + sR * (uMR - uR));
            }
            else if (sR < 0.0)
                return fR;
        }

    private:
        std::shared_ptr<Model> model;
};

#endif // HYPSYS1D_NUMERICAL_FLUX_HPP
