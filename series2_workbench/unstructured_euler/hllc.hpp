#pragma once

// calculating left and right wave speeds sL, sR according to Einfeldt (1988):
std::pair<double, double> calculate_sL_sR(const EulerState& uL, const EulerState& uR)
{
    // initial implementation with simple speed estimates (NOT Einfeldt):
    return {euler::minEigenValue(uL, uR), euler::maxEigenValue(uL, uR)};
}

// calculate intermediate wave speed s_star:
double s_star(const EulerState& uL, const EulerState& uR, double sL, double sR)
{
    auto pL= euler::pressure(uL);
    auto pR= euler::pressure(uR);

    auto rhoL= uL[0];
    auto rhoR= uR[0];

    auto vL= std::sqrt(uL[1] * uL[1] + uL[2] * uL[2]) / rhoL;
    auto vR= std::sqrt(uR[1] * uR[1] + uR[2] * uR[2]) / rhoR;

    return (pR - pL + rhoL * vL * (sL - vL) - rhoR * vR * (sR - vR))
           / (rhoL * (sL - vL) - rhoR * (sR - vR));
}

// calculate intermediate left/right state u_star:
EulerState u_star(const EulerState& uK, double sK, double s_star)
{
    auto rho= uK[0];
    auto m1= uK[1];
    auto m2= uK[2];
    auto E= uK[3];
    auto u1= m1 / rho;
    auto u2= m2 / rho;
    auto p= euler::pressure(uK);

    EulerState uM;

    uM[0]= 1.0;
    uM[1]= s_star;
    uM[2]= u2;
    uM[3]= E / rho + (s_star - u1) * (s_star + p / (rho * (sK - u1)));

    return rho * (sK - u1) / (sK - s_star) * uM;
}

// calculate intermediate left/right flux flux_star:
EulerState flux_star(const EulerState& u_star, const EulerState& uK, double sK)
{
    return euler::flux(uK) + sK * (u_star - uK);
}

/// HLLC numerical flux with Einfeldt-Batten wavespeeds.
/** Reference: Batten, Wavespeed Estimates for the HLLC Riemann Solver, 1997
 *  @param euler
 *  @param uL    conserved variables 'inside' of the cell.
 *  @param uR    conserved variables 'outside' of the cell.
 */
EulerState hllc(const EulerState &uL, const EulerState &uR) {
    // implement the HLLC approximate Riemann solver.
    // tip, compute the three wave speeds sL, s_star & sR in
    // a separate function to keep things readable.

    auto [sL, sR]= calculate_sL_sR(uL, uR);

    if (sL > 0)
        return euler::flux(uL);
    
    if (sR < 0)
        return euler::flux(uR);

    auto sM= s_star(uL, uR, sL, sR);

    if (sM > 0) // K = L: left intermediate flux
        return flux_star(u_star(uL, sL, sM), uL, sL);

    // if (sM <= 0) implicitly: K = R: right intermediate flux
    return flux_star(u_star(uR, sR, sM), uR, sR);
}
