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
                       H(u) - (v(u) * c(u)),   0.5 * v(u) * v(u),   H(u) - (v(u) * c(u));

            return eigvecs;
        }
       
        double max_eigenvalue(const Eigen::VectorXd &u) const override
        {
            // check if we want cwiseAbs here:
            return (eigenvalues(u).cwiseAbs()).maxCoeff();
            // return (eigenvalues(u).maxCoeff());
        }

        Eigen::VectorXd cons_to_prim(const Eigen::VectorXd &u_cons) const override
        {
            Eigen::VectorXd u_prim(n_vars);
            
            double m, p, rho, v, E;
            rho= u_cons(0);
            m= u_cons(1);
            E= u_cons(2);
            v= m / rho;
            p= (gamma - 1) * (E - 0.5 * rho * v * v);

            u_prim << rho,
                      v,
                      p;
                      
            return u_prim;
        }
        
        Eigen::VectorXd prim_to_cons(const Eigen::VectorXd &u_prim) const override
        {
            Eigen::VectorXd u_cons(n_vars);

            u_cons << rho(u_prim),
                      m(u_prim),
                      E(u_prim);

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
        /// Vector u_prim must contain primary variables: u_prim = (rho, v, p)
        /// ------------------------------------------------------------------

        /// Primitive variables:
        // rho : density
        inline double rho(const Eigen::VectorXd& u_prim) const override
        {
            return u_prim(0);
        }

        // v : velocity
        inline double v(const Eigen::VectorXd& u_prim) const override
        {
            return u_prim(1);
        }

        // p : pressure
        inline double p(const Eigen::VectorXd& u_prim) const override
        {
            return u_prim(2);
        }



        /// Conserved variables:
        // m : momentum
        inline double m(const Eigen::VectorXd& u_prim) const override
        {
            return rho(u_prim) * v(u_prim);
        }

        // E : total energy for ideal polytropic gas (internal energy + kinetic energy)
        inline double E(const Eigen::VectorXd& u_prim) const override
        {
            return p(u_prim) / (gamma - 1) + 0.5 * rho(u_prim) * v(u_prim) * v(u_prim);
        }


        /// c : speed of sound
        inline double c(const Eigen::VectorXd& u_prim) const override
        {
            return std::sqrt(gamma * p(u_prim) / rho(u_prim));
        }

        /// H : total specific enthalpy
        inline double H(const Eigen::VectorXd& u_prim) const override
        {
            return (E(u_prim) + p(u_prim)) / rho(u_prim);
        }
        
};


std::shared_ptr<Model> make_model (const nlohmann::json &config);

#endif // HYPSYS1D_MODEL_HPP
