#include "../../include/ancse/polynomial_basis.hpp"

/// Computes the Legendre polynomial basis
/// at a given reference point xi \in [0,1]
Eigen::VectorXd PolynomialBasis::operator() (double xi) const
{
    // coefficient vector for degree p:
    Eigen::VectorXd phi(p + 1);

    // lambda to get first three Legendre polynomials:
    auto legendre= [](int k, double x) -> double
    {
        if (k == 0)
            return 1.0;
        else if (k == 1)
            return x;
        else if (k == 2)
            return 0.5 * (3.0 * x * x - 1.0);
        else
            std::abort();
    };

    for (int k= 0; k < p + 1; ++k)
        phi(k)= std::sqrt(2 * k + 1) * legendre(k, 2 * xi - 1);

    return scaling_factor * phi;
}

/// Computes the derivative of Legendre polynomial basis
/// at a given reference point xi \in [0,1]
Eigen::VectorXd PolynomialBasis::deriv (double xi) const
{
    // spatial derivative coefficient vector for degree p:
    Eigen::VectorXd derivPhi(p + 1);

    // lambda to get first three Legendre polynomial derivatives:
    auto deriv_legendre= [](int k, double x) -> double
    {
        if (k == 0)
            return 0.0;
        else if (k == 1)
            return 1.0;
        else if (k == 2)
            return 3.0 * x;
        else
            std::abort();
    };

    for (int k= 0; k < p + 1; ++k)
        derivPhi(k)= std::sqrt(2 * k + 1) * deriv_legendre(k, 2 * xi - 1);
    
    return scaling_factor * derivPhi;
}
