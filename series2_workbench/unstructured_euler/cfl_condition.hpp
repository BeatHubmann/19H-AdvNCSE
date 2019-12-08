#pragma once

#include <fmt/format.h>

#include "euler.hpp"
#include "mesh.hpp"

class CFLCondition {
  public:
    explicit CFLCondition(const Mesh &mesh) : dx(mesh.getMinimumInradius()) {}

    double operator()(const Eigen::MatrixXd &U) const {
        // compute the cfl condition here.
        // you can use `assert_valid_timestep` to check if
        // the computed value is valid.
        
        // in analogy to exercise 1 and according to eqns 33, 34 of exercise sheet:
        
        const int n_cells= U.cols();

        double a_max= 0.0;

        for (int i= 0; i < n_cells; ++i)
            a_max= std::max(a_max, euler::maxEigenValue(U.col(i)));

        double dt= cfl_number * dx / a_max;
        assert_valid_timestep(dt);

        return dt;
    }

    void assert_valid_timestep(double dt_cfl) const {
        if (dt_cfl <= 0.0 || !std::isfinite(dt_cfl)) {
            throw std::runtime_error(
                fmt::format("Non-positive timestep: dt = {:.3e}", dt_cfl));
        }
    }

  private:
    double dx;
    double cfl_number = 0.45;
};