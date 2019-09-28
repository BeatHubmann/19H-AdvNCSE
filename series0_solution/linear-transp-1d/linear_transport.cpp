#include "writer.hpp"
#include <Eigen/Core>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>

void apply_boundary_conditions(Eigen::ArrayXXd &u, int k) {
    auto N = u.rows();
    u(0, k) = u(1, k);
    u(N - 1, k) = u(N - 2, k);
}

//----------------upwindFDBegin----------------
/// Uses forward Euler and upwind finite differences to compute u from time 0 to
/// time T
///
/// @param[in] u0 the initial conditions in the physical domain, excluding
/// ghost-points.
/// @param[in] dt the time step size
/// @param[in] T the solution is computed for the interval [0, T]. T is assumed
/// to be a multiple of of dt.
/// @param[in] a the advection velocity
/// @param[in] domain left & right limit of the domain
///
/// @return returns the solution 'u' at every time-step and the corresponding
/// time-steps. The solution `u` includes the ghost-points.
std::pair<Eigen::ArrayXXd, Eigen::VectorXd>
upwindFD(const Eigen::VectorXd &u0,
         double dt,
         double T,
         const std::function<double(double)> &a,
         const std::pair<double, double> &domain) {

    auto N = u0.size();
    auto nsteps = int(round(T / dt));

    auto u = Eigen::ArrayXXd(N + 2, nsteps + 1);
    auto time = Eigen::VectorXd(nsteps + 1);

    auto [xL, xR] = domain;
    double dx = (xR - xL) / (N - 1.0);

    /* Initialize u */
    //// ANCSE_START_TEMPLATE
    u.col(0).segment(1, N) = u0;
    apply_boundary_conditions(u, 0);
    time[0] = 0.0;
    //// ANCSE_END_TEMPLATE

    /* Main loop */
    //// ANCSE_START_TEMPLATE
    for (int k = 0; k < nsteps; k++) {
        for (int j = 1; j < N + 1; j++) {

            double x = xL + (j - 1) * dx;

            double uL = u(j - 1, k);
            double uM = u(j, k);
            double uR = u(j + 1, k);
            double ax = a(x);
            double c = 0.5 * dt / dx;

            u(j, k + 1)
                = uM - c * (ax * (uR - uL) - fabs(ax) * (uR - 2 * uM + uL));
        }
        /* Outflow boundary conditions */
        apply_boundary_conditions(u, k + 1);
        time[k + 1] = (k + 1) * dt;
    }
    //// ANCSE_END_TEMPLATE

    return {std::move(u), std::move(time)};
}
//----------------upwindFDEnd----------------

//----------------centeredFDBegin----------------
/// Uses forward Euler and centered finite differences to compute u from time 0
/// to time T
///
/// @param[in] u0 the initial conditions, as column vector
/// @param[in] dt the time step size
/// @param[in] T the solution is computed for the interval [0, T]. T is assumed
/// to be a multiple of of dt.
/// @param[in] a the advection velocity
/// @param[in] domain left & right limit of the domain
///
/// @return returns the solution 'u' at every time-step and the corresponding
/// time-steps. The solution `u` includes the ghost-points.
std::pair<Eigen::ArrayXXd, Eigen::VectorXd>
centeredFD(const Eigen::VectorXd &u0,
           double dt,
           double T,
           const std::function<double(double)> &a,
           const std::pair<double, double> &domain) {

    auto N = u0.size();
    auto nsteps = int(round(T / dt));
    auto u = Eigen::ArrayXXd(N + 2, nsteps + 1);
    auto time = Eigen::VectorXd(nsteps + 1);

    auto [xL, xR] = domain;
    double dx = (xR - xL) / (N - 1.0);

    /* Initialize u */
    //// ANCSE_START_TEMPLATE
    u.col(0).segment(1, N) = u0;
    apply_boundary_conditions(u, 0);
    time[0] = 0.0;
    //// ANCSE_END_TEMPLATE

    /* Main loop */
    //// ANCSE_START_TEMPLATE
    for (int k = 0; k < nsteps; k++) {
        for (int j = 1; j < N + 1; j++) {
            double x = xL + (j - 1) * dx;
            u(j, k + 1)
                = u(j, k)
                  - dt / (2.0 * dx) * a(x) * (u(j + 1, k) - u(j - 1, k));
        }

        /* Outflow boundary conditions */
        apply_boundary_conditions(u, k + 1);
        time[k + 1] = (k + 1) * dt;
    }
    //// ANCSE_END_TEMPLATE

    return {std::move(u), std::move(time)};
}
//----------------centeredFDEnd----------------

/* Initial condition: rectangle */
double ic(double x) {
    if (x < 0.25 || x > 0.75)
        return 0.0;
    else
        return 2.0;
}

int main() {
    double T = 2.0;
    double dt = 0.002; // Change this for timestep comparison
    int N = 101;

    double xL = 0.0;
    double xR = 5.0;
    auto domain = std::pair<double, double>{xL, xR};

    auto a = [](double x) { return std::sin(2.0 * M_PI * x); };

    Eigen::VectorXd u0(N);
    double h = (xR - xL) / (N - 1.0);
    /* Initialize u0 */
    for (int i = 0; i < u0.size(); i++) {
        u0[i] = ic(xL + h * i);
    }

    const auto &[u_upwind, time_upwind] = upwindFD(u0, dt, T, a, domain);
    writeToFile("time_upwind.txt", time_upwind);
    writeMatrixToFile("u_upwind.txt", u_upwind);

    const auto &[u_centered, time_centered] = centeredFD(u0, dt, T, a, domain);
    writeToFile("time_centered.txt", time_centered);
    writeMatrixToFile("u_centered.txt", u_centered);
}
