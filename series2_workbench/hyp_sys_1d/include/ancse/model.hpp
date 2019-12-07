#ifndef HYPSYS1D_MODEL_HPP
#define HYPSYS1D_MODEL_HPP

#include <cmath>
#include <memory>

#include <Eigen/Dense>
#include <ancse/config.hpp>

/// Interface for implementing different models,
/// eg. Euler equations, Shallow-water equations
///
/// Add more functions to this interface if needed.
class Model {
  public:
    virtual ~Model() = default;

    // Pure virtual functions:
    virtual Eigen::VectorXd flux(const Eigen::VectorXd& u) const = 0;
    virtual Eigen::VectorXd eigenvalues(const Eigen::VectorXd& u) const = 0;
    virtual Eigen::MatrixXd eigenvectors(const Eigen::VectorXd& u) const = 0;
    virtual double max_eigenvalue(const Eigen::VectorXd& u) const = 0;

    virtual Eigen::VectorXd cons_to_prim(const Eigen::VectorXd& u) const = 0;
    virtual Eigen::VectorXd prim_to_cons(const Eigen::VectorXd& u) const = 0;

    virtual int get_nvars() const = 0;
    virtual std::string get_name() const = 0;

    // Virtual functions:
    virtual std::pair<double, double> lo_hi_eigenvalues(const Eigen::VectorXd& u) const {}

    virtual double rho(const Eigen::VectorXd& u_prim) const {}
    virtual double v(const Eigen::VectorXd& u_prim) const {}
    virtual double p(const Eigen::VectorXd& u_prim) const {}
    virtual double m(const Eigen::VectorXd& u_prim) const {}
    virtual double E(const Eigen::VectorXd& u_prim) const {}
    virtual double c(const Eigen::VectorXd& u_prim) const {}
    virtual double H(const Eigen::VectorXd& u_prim) const {}

    virtual void set_gamma(const double gamma_) {}
    virtual double get_gamma() const {}
};

class Burgers : public Model {
  public:

    Eigen::VectorXd flux(const Eigen::VectorXd &u) const override
    {
        Eigen::VectorXd f(n_vars);
        f(0) = 0.5*u(0)*u(0);

        return f;
    }

    Eigen::VectorXd eigenvalues(const Eigen::VectorXd &u) const override
    {
        Eigen::VectorXd eigvals(n_vars);
        eigvals(0) = u(0);

        return eigvals;
    }

    Eigen::MatrixXd eigenvectors(const Eigen::VectorXd &) const override
    {
        Eigen::MatrixXd eigvecs(n_vars, n_vars);
        eigvecs (0,0) = 1;

        return eigvecs;
    }

    double max_eigenvalue(const Eigen::VectorXd &u) const override {
        return (eigenvalues(u).cwiseAbs()).maxCoeff();
    }

    Eigen::VectorXd cons_to_prim(const Eigen::VectorXd &u) const override {
        return u;
    }

    Eigen::VectorXd prim_to_cons(const Eigen::VectorXd &u) const override {
        return u;
    }

    int get_nvars() const override
    {
        return n_vars;
    }

    std::string get_name() const override
    {
        return "burgers";
    }

  private:
    int n_vars = 1;
};

/// Euler equations
class Euler : public Model {
    public:

        Eigen::VectorXd flux(const Eigen::VectorXd &u) const override
        {
            Eigen::VectorXd f(n_vars);
            
            f << rho(u) * v(u),
                 rho(u) * v(u) * v(u) + p(u),
                 (E(u) + p(u)) * v(u);
                
            return f;
        }
        
        Eigen::VectorXd eigenvalues(const Eigen::VectorXd &u) const override
        {
            Eigen::VectorXd eigvals(n_vars);
        
            eigvals << v(u) - c(u),
                       v(u),
                       v(u) + c(u);

            return eigvals;
        }

        Eigen::MatrixXd eigenvectors(const Eigen::VectorXd &u) const override
        {
            Eigen::MatrixXd eigvecs(n_vars, n_vars);
            
            eigvecs << 1.0,                    1.0,                 1.0,
                       v(u) - c(u),            v(u),                v(u) + c(u),
                       H(u) - v(u) * c(u),     0.5 * v(u) * v(u),   H(u) + (v(u) * c(u));

            return eigvecs;
        }
       
        double max_eigenvalue(const Eigen::VectorXd &u) const override
        {
            // actually max(abs(eigenvalues)):
            return (eigenvalues(u).cwiseAbs()).maxCoeff();
        }

        Eigen::VectorXd cons_to_prim(const Eigen::VectorXd &u_cons) const override
        {
            Eigen::VectorXd u_prim(n_vars);
            
            u_prim << rho(u_cons),
                      v(u_cons),
                      p(u_cons);
                      
            return u_prim;
        }
        
        Eigen::VectorXd prim_to_cons(const Eigen::VectorXd &u_prim) const override
        {
            Eigen::VectorXd u_cons(n_vars);

            u_cons << u_prim(0),
                      u_prim(0) * u_prim(1),
                      u_prim(2) / (gamma - 1) + 0.5 * u_prim(0) * u_prim(1) * u_prim(1);

            return u_cons;
        }

        inline void set_gamma(const double gamma_) override
        {
            gamma= gamma_;
        }

        inline double get_gamma() const override
        {
            return gamma;
        }

        inline int get_nvars() const override
        {
            return n_vars;
        }

        std::string get_name() const override
        {
            return name;
        }

    private:
        // 3D:
        const int n_vars= 3;
        // Monatomic gas:
        double gamma= 5./3.;

        inline static const std::string name= "euler";


        std::pair<double, double> lo_hi_eigenvalues(const Eigen::VectorXd &u) const override
        {
            // check if we want cwiseAbs here:
            // return (eigenvalues(u).cwiseAbs()).minCoeff();
            Eigen::VectorXd eigenvals= eigenvalues(u);
            return {eigenvals.minCoeff(), eigenvals.maxCoeff()};
        }

        /// Helper functions to deal w/ common names for Euler eqn expressions
        /// Vector u must contain conserved variables: u = (rho, m, E)
        /// as given in tasks initial conditions
        /// ------------------------------------------------------------------

        /// Conserved variables:
        // rho : density
        inline double rho(const Eigen::VectorXd& u) const override
        {
            return u(0);
        }

        // m : momentum
        inline double m(const Eigen::VectorXd& u) const override
        {
            return u(1);
        }
        
        // E : total energy for ideal polytropic gas (internal energy + kinetic energy)
        inline double E(const Eigen::VectorXd& u) const override
        {
            return u(2);
        }

        /// Primitive variables:
        // rho : density
        // same as for conserved variables

        // v : velocity
        inline double v(const Eigen::VectorXd& u) const override
        {
            return m(u) / rho(u);
        }

        // p : pressure
        inline double p(const Eigen::VectorXd& u) const override
        {
            return (E(u) - 0.5 * m(u) * m(u) / rho(u)) * (gamma - 1);
        }

        /// c : speed of sound
        inline double c(const Eigen::VectorXd& u) const override
        {
            return std::sqrt(gamma * p(u) / rho(u));
        }

        /// H : total specific enthalpy
        inline double H(const Eigen::VectorXd& u) const override
        {
            return (E(u) + p(u)) / rho(u);
        }
};

std::shared_ptr<Model> make_model (const nlohmann::json &config);

#endif // HYPSYS1D_MODEL_HPP