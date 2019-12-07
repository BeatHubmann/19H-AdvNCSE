#ifndef HYPSYS1D_LIMITERS_HPP
#define HYPSYS1D_LIMITERS_HPP

#include <cmath>
#include <algorithm>

// Slope limiter helpers:
template <typename T>
inline T slope(T left, T right, T h)
{
    return (right - left) / h;
}

template <typename T>
inline T sign(T a)
{
    return copysign((T)1.0, a);
}

template <typename T>
inline T minmod(T a, T b)
{
    return 0.5 * (sign(a) + sign(b)) * std::min(std::abs(a), std::abs(b));
}

template <typename T>
inline T maxmod(T a, T b)
{
    return 0.5 * (sign(a) + sign(b)) * std::max(std::abs(a), std::abs(b));
}

template <typename T>
inline T minmod(T a, T b, T c)
{
    return minmod(a, minmod(b, c));
}

template <typename T>
inline T minabs(T a, T b)
{
    const double tolerance= 1e-10;
    if ((double) (std::abs(std::abs(a) - std::abs(b))) < tolerance)
        return 0.5*(a+b);
    return std::abs(a) < std::abs(b) ? a : b;
}

/// FVM slope limiters
struct MinMod
{
    inline double operator()(double sL, double sR) const
    {
        return minmod(sL, sR);
    }
};

struct SuperBee
{
    inline double operator()(double sL, double sR) const
    {
        double A = minmod(2.0 * sL, sR);
        double B = minmod(sL, 2.0 * sR);

        return maxmod(A, B);
    }
};

struct MonotonizedCentral
{
    inline double operator()(double sL, double sR) const
    {
        return minmod(2.0 * sL, 0.5 * (sL + sR), 2.0 * sR);
    }
};

struct MinAbs
{
    inline double operator()(double sL, double sR) const
    {
        return minabs(sL, sR);
    }
};


struct VanLeerFVM
{
    inline double operator()(double sL, double sR) const
    {
        double r = sL / (sR + eps);
        return (r + std::abs(r)) / (1 + std::abs(r)) * sR;
    }

    private:
        double eps = 1e-10;
};


struct Unlimited
{
    inline double operator()(double a, double b) const
    {
        return 0.5 * (a + b);
    }
};

/// DG limiters
struct VanLeer
{
    inline double operator()(double s, double sm, double sp) const
    {
        return minmod(s, sm, sp);
    }
};

struct Shu
{
    explicit Shu (const double dx_) : dx (dx_) {}

    inline double operator()(double s, double sm, double sp) const
    {
        if (std::abs(s) < M * dx * dx)
        // {
            return s;
        // }
        // else
        // {
        return minmod(s, sm, sp);
        // }
    }

    private:
        double dx;
        double M = 50;
};


#endif // HYPSYS1D_LIMITERS_HPP
