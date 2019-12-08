#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

#include "gradient.hpp"
#include "hllc.hpp"
#include "mesh.hpp"
#include "slope_limiter.hpp"

// Note: this class will compute the rate of change due to the fluxes.
// Note: the reason we made this a class is that it allows you to allocate
//       buffers, once at the beginning of the simulation. Add these buffers
//       as needed.
class FluxRateOfChange {
  public:
    explicit FluxRateOfChange(int n_cells) : n_cells(n_cells) {}

    void operator()(Eigen::MatrixXd &dudt,
                    const Eigen::MatrixXd &u,
                    const Mesh &mesh) const
    {
        // Compute the rate of change of u.
        // Note: Please use the method `computeFlux` to abstract
        // away the details of computing the flux through a
        // given interface.
        // Note: You can use `assert_valid_flux` to check
        // if what `computeFlux` returns makes any sense.
        // Note: Do not assume `dudt` is filled with zeros.

        // sanity check:
        assert ((mesh.getNumberOfTriangles() == n_cells));
        
        // reset dudt to zero:
        dudt*= 0.0;

        // iterate over cells then over the 3 edges within each cell:
        for (int i= 0; i < n_cells; ++i)
            for (int k= 0; k < 3; ++k)
            {
                EulerState edge_flux= computeFlux(u, i, k, mesh);
                assert_valid_flux(mesh, i, k, edge_flux);
                dudt.row(i)+= edge_flux;
            }
    }

    void assert_valid_flux(const Mesh &mesh,
                           int i,
                           int k,
                           const EulerState &nF) const {
        // This is mostly for debugging (but also important to check in
        // real simulations!): Make sure our flux contribution is not
        // nan (ie. it is not not a number, ie it is a number)
        if (!euler::isValidFlux(nF)) {
            // clang-format off
            throw std::runtime_error(
                "invalid value detected in numerical flux, " + euler::to_string(nF)
                + "\nat triangle: " + std::to_string(i)
                + "\nedge:        " + std::to_string(k)
                + "\nis_boundary: " + std::to_string(!mesh.isValidNeighbour(i, k)));
            // clang-format on
        }
    }

    /// Compute the flux through the k-th interface of cell i.
    EulerState computeFlux(const Eigen::MatrixXd &U,
                           int i,
                           int k,
                           const Mesh &mesh) const
    {
        auto boundary_type = mesh.getBoundaryType(i, k);

        if (boundary_type == Mesh::BoundaryType::INTERIOR_EDGE)
        {
            return computeInteriorFlux(U, i, k, mesh);
        }
        else
        {
            if (boundary_type == Mesh::BoundaryType::OUTFLOW_EDGE)
            {
                return computeOutflowFlux(U, i, k, mesh);
            }
            else /* boundary_type == Mesh::BoundaryType::WING_EDGE */
            {
                return computeReflectiveFlux(U, i, k, mesh);
            }
        }
    }

    /// Compute the outflow flux through the k-th interface of cell i.
    /** Note: you know that edge k is an outflow edge.
     */
    EulerState computeOutflowFlux(const Eigen::MatrixXd &U,
                                  int i,
                                  int k,
                                  const Mesh &mesh) const
    {
        // Implement the outflow flux boundary condition.
        auto u_aligned= rotate_align(U.row(i), i, k, mesh);

        auto f_aligned= euler::flux(u_aligned);

        auto f= rotate_dealign(f_aligned, i, k, mesh);

        return mesh.getEdgeLength(i, k) * f; // eqn (35)
    }

    /// Compute the reflective boundary flux through the k-th edge of cell i.
    /** Note: you know that edge k is a reflective/wall boundary edge.
     */
    EulerState computeReflectiveFlux(const Eigen::MatrixXd &U,
                                     int i,
                                     int k,
                                     const Mesh &mesh) const
    {
        // Implement the reflective flux boundary condition.
        auto u= U.row(i);

        // get unit outward normal and transverse t:
        auto n= mesh.getUnitNormal(i, k).normalized();
        auto t= Eigen::Vector2d(-n(1), n(0));

        // assemble u_star:
        double rho= u[0];
        Eigen::Vector2d v= u.segment(1, 2);
        double E= u[3];

        EulerState u_star;
        u_star[0]=  rho;
        u_star.segment(1, 2)= -rho * v.dot(n) * n + rho * v.dot(t) * t;
        u_star[3]=  E;

        // rotate:
        auto u_aligned= rotate_align(u, i, k, mesh);
        auto u_star_aligned= - rotate_align(u_star, i, k, mesh); // (-) for reflection
        
        // flux:
        auto f_aligned= hllc(u_aligned, u_star_aligned); // equation (36)

        // derotate:
        auto f= rotate_dealign(f_aligned, i, k, mesh);
        
        return mesh.getEdgeLength(i, k) * f;
    }

    /// Compute the flux through the k-th interface of cell i.
    /** Note: This edge is an interior edge, therefore approximate the flux
     * through this edge with the appropriate FVM formulas.
     */
    EulerState computeInteriorFlux(const Eigen::MatrixXd &U,
                                   int i,
                                   int k,
                                   const Mesh &mesh) const
    {
        // Reconstruct the trace values of U and compute
        // the numerical flux through the k-th interface of
        // cell i.

        // figure out neighbour of i at edge k and check validity:
        int j= mesh.getNeighbour(i, k);
        assert(j >= 0);

        // figure out neighbour j's edge l collocated with i's edge k:
        int l= mesh.getNeighbourEdge(i, k);

        // reconstruction:
        auto uL= reconstruction(U.row(i), i, k, mesh);
        auto uR= reconstruction(U.row(j), j, l, mesh);

        // rotate:
        auto uL_aligned=   rotate_align(uL, i, k, mesh);
        auto uR_aligned= - rotate_align(uR, j, l, mesh);

        // flux:
        auto f_aligned= hllc(uL_aligned, uR_aligned); // equation (32)

        // derotate:
        auto f= rotate_dealign(f_aligned, i, k, mesh);
        
        return mesh.getEdgeLength(i, k) * f;
    }


    private:
        // generate rotation matrix to align edge k of triangle i with cartesian y-axis:
        Eigen::Matrix2d rot(int i, int k, const Mesh& mesh) const
        {
            auto n= mesh.getUnitNormal(i, k).normalized(); // why is this called UnitNormal if it's not unit normalized?

            Eigen::Matrix2d rot;
            rot << n(0),  n(1),
                  -n(1),  n(0);

            return rot;
        }

        // generate inverse rotation matrix to DEalign edge k of triangle i from cartesian y-axis back to original position:
        Eigen::Matrix2d rot_inv(int i, int k, const Mesh& mesh) const
        {
            // our rotation matrix is orthogonal: RR'=I:
            return rot(i, k, mesh).transpose();
        }

        // rotate u to align velocity component [u(1), u(2)] with cartesian x-axis:
        EulerState rotate_align(EulerState u, int i, int k, const Mesh& mesh) const
        {
            u.segment(1, 2)= rot(i, k, mesh) * u.segment(1, 2);
            return u;
        }

        // rotate calculated f back to point in original grid position before rotation:
        EulerState rotate_dealign(EulerState f, int i, int k, const Mesh& mesh) const
        {
            f.segment(1, 2)= rot_inv(i, k, mesh) * f.segment(1, 2);
        }

        // arggh there's probably no time for doing piecewise linear construction in task 2h)
        EulerState reconstruction(const Eigen::MatrixXd& U, int i, int k, const Mesh& mesh) const
        {
            // pass:
            return U.row(i);
        }

        int n_cells;
};
