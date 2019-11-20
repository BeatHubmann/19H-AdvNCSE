#ifndef FVMSCALAR1D_RECONSTRUCTION_HPP
#define FVMSCALAR1D_RECONSTRUCTION_HPP

#include <Eigen/Dense>
#include <cmath>
#include <memory>

#include <ancse/grid.hpp>
#include <ancse/model.hpp>
#include <ancse/rate_of_change.hpp>
#include <ancse/simulation_time.hpp>

inline double sign(double a) { return copysign(1.0, a); }

//// ANCSE_CUT_START_TEMPLATE
//----------------SlopeLimiterABegin----------------
inline double minmod(double a, double b) {
    return 0.5 * (sign(a) + sign(b)) * std::min(std::abs(a), std::abs(b));
}

inline double maxmod(double a, double b) {
    return 0.5 * (sign(a) + sign(b)) * std::max(std::abs(a), std::abs(b));
}

inline double minmod(double a, double b, double c) {
    return minmod(a, minmod(b, c));
}

inline double minabs(double a, double b) {
    double tolerance = 1e-10;
    if (std::abs( std::abs(a) - std::abs(b) ) < tolerance)
        return 0.5*(a+b);
    return std::abs(a) < std::abs(b) ? a : b;
}
//----------------SlopeLimiterAEnd----------------

//----------------SlopeLimiterBBegin----------------
struct MinMod {
    inline double operator()(double sL, double sR) const {
        return minmod(sL, sR);
    }
};
//----------------SlopeLimiterBEnd----------------

//----------------SlopeLimiterCBegin----------------
struct MinAbs {
    inline double operator()(double sL, double sR) const {
        return minabs(sL, sR);
    }
};

struct SuperBee {
    inline double operator()(double sL, double sR) const {
        double A = minmod(2.0 * sL, sR);
        double B = minmod(sL, 2.0 * sR);

        return maxmod(A, B);
    }
};

struct VanLeer {
    inline double operator()(double sL, double sR) const {
        double r = sL / (sR + eps);
        return (r + std::abs(r)) / (1 + std::abs(r)) * sR;
    }

  private:
    double eps = 1e-10;
};

struct MonotonizedCentral {
    inline double operator()(double sL, double sR) const {
        return minmod(2.0 * sL, 0.5 * (sL + sR), 2.0 * sR);
    }
};

struct Unlimited {
    inline double operator()(double a, double b) const { return 0.5 * (a + b); }
};
//----------------SlopeLimiterCEnd----------------
//// ANCSE_END_TEMPLATE

class PWConstantReconstruction {
  public:
    /// Compute the left and right trace at the interface i + 1/2.
    /** Note: This API is agnostic to the number of cell-averages required
     *        by the method. Therefore, reconstructions with different stencil
     *        sizes can implement this API; and this call can be used in parts
     *        of the code that do not need to know about the details of the
     *        reconstruction.
     */
    inline std::pair<double, double> operator()(const Eigen::VectorXd &u,
                                                int i) const {
        return (*this)(u[i], u[i + 1]);
    }

    /// Compute the left and right trace at the interface.
    /** Piecewise constant reconstruction of the left and right trace only
     *  requires the cell-average to the left and right of the interface.
     *
     *  Note: Compared to the other overload this reduces the assumption on
     *        how the cell-averages are stored. This is useful when testing and
     *        generally makes the function useful in more situations.
     */
    inline std::pair<double, double> operator()(double ua, double ub) const {
        return {ua, ub};
    }
};

//// ANCSE_CUT_START_TEMPLATE
//----------------LinearRCBegin----------------
template <class SlopeLimiter>
class PWLinearReconstruction {
  public:
    explicit PWLinearReconstruction(const SlopeLimiter &slope_limiter)
        : slope_limiter(slope_limiter) {}

    std::pair<double, double> operator()(const Eigen::VectorXd &u,
                                         int i) const {
        return (*this)(u[i - 1], u[i], u[i + 1], u[i + 2]);
    }

    std::pair<double, double>
    operator()(double ua, double ub, double uc, double ud) const {
        double sL = ub - ua;
        double sM = uc - ub;
        double sR = ud - uc;

        double uL = ub + 0.5 * slope_limiter(sL, sM);
        double uR = uc - 0.5 * slope_limiter(sM, sR);

        return {uL, uR};
    }

  private:
    SlopeLimiter slope_limiter;
};
//----------------LinearRCEnd----------------
//// ANCSE_END_TEMPLATE

#endif // FVMSCALAR1D_RATE_OF_CHANGE_HPP
