#include "writer.hpp"
#include <Eigen/Core>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <nlohmann/json.hpp>
#include <filesystem>

//! Initial condition: grid-aligned rectangle
//! @param[in] (x,y) point to evaulate the initial condition
//! @param[out] evaluation of the function at (x,y)
double ic(double x, double y) {
    if (x > -0.3 && x < 0.3 && y > -0.3 && y < 0.3)
        return 1.0;
    else
        return 0.0;
}

//! Load configuration values from a JSON file
//! @param[in] jsonFilename file (including path if needed) of config.json
//! @param[out] json object with parameters: T, {x,y}{min,max}, Nx, Ny, cfl
nlohmann::json loadConfig(std::string jsonFilename) {
    std::ifstream i(jsonFilename);
    assert(i.good() && "config.json not found in current or parent directory");

    nlohmann::json j;
    i >> j;
    return j;
}

//! Apply periodic boundary conditions to matrix u.
//! u has relevant values in the Nx x Ny submatrix at the center
//! @param[in] u  (Nx+2)x(Ny+2) matrix
void applyBoundaryConditions(Eigen::MatrixXd &u) {
    //// ANCSE_START_TEMPLATE
    //----------------BCStart----------------
    Eigen::Index r = u.rows();
    Eigen::Index c = u.cols();

    for(int i = 1 ; i < r-1 ; i++) {
        u(i, c-1) = u(i, 1);
        u(i, 0) = u(i, c-2);
    }
    for( int j = 1 ; j < c-1; j++) {
        u(r-1, j) = u(1, j);
        u(0, j) = u(r-2, j);
    }
    // this is first order, so we ignore corners
    //----------------BCEnd----------------
    //// ANCSE_END_TEMPLATE
}

//! An implementation of the 1D upwind numerical flux
//! @param[in] uM  value to the left of the interface
//! @param[in] uR  value to the right of the interface
//! @param[in] a   velocity at the interface
double F(double uM, double uP, double a) {
    return std::max(a, 0.)*uM + std::min(a, 0.)*uP;
}

//! Compute one step of the upwind method in 2D
//! @param[in] u           matrix of size (Nx+2)x(Ny+2) which will contain U^{n+1}
//! @param[in] u_old       matrix of size (Nx+2)x(Ny+2) of values U^{n}
//! @param[in] dx, dy, dt  meshsteps and timestep
//! @param[in] a           velocity as a function of R^2 to R^2
void updateUpwind(Eigen::MatrixXd &u, Eigen::MatrixXd &u_old, double dx,
                  double dy, double dt, double xmin, double ymin,
                  const std::function<Eigen::Vector2d(double, double)> &a) {
    //// ANCSE_START_TEMPLATE
    //----------------UpdateStart----------------
    for (int i = 1 ; i < u_old.cols()-1 ; i++) {
        for (int j = 1 ; j < u_old.rows()-1 ; j++) {
            auto x_ctr = xmin + (i-0.5)*dx;
            auto y_ctr = ymin + (j-0.5)*dy;
            auto a_N = a(x_ctr, y_ctr + dy/2);
            auto a_S = a(x_ctr, y_ctr - dy/2);
            auto a_E = a(x_ctr + dx/2, y_ctr);
            auto a_W = a(x_ctr - dx/2, y_ctr);
            u(i,j) = u_old(i,j) - (dt/dy)*(F(u_old(i,j), u_old(i, j+1), a_N(1)) - F(u_old(i,j-1), u_old(i, j), a_S(1)))
                                - (dt/dx)*(F(u_old(i,j), u_old(i+1, j), a_E(0)) - F(u_old(i-1,j), u_old(i, j), a_W(0))) ;
        }
    }
    //----------------UpdateEnd----------------
    //// ANCSE_END_TEMPLATE
}

//! Clear the contents of the output file before starting
//! Useful because we write in append mode
//! @param[in] outfile  name of the file to be wiped
void wipeFile(std::string outfile) {
    std::ofstream outstrm;
    outstrm.open(outfile, std::ofstream::out | std::ofstream::trunc);
    outstrm.close();
}

int main() {
    // Path to the config file relative to the binary. Try a couple of likely locations.
    std::string config_file = std::filesystem::exists("../config.json") ? "../config.json" : "config.json";

    auto j = loadConfig(config_file);
    double T = j["T"];
    double xmin = j["xmin"], ymin = j["ymin"], xmax = j["xmax"], ymax = j["ymax"];
    int Nx = j["Nx"], Ny = j["Ny"];
    double cfl = j["cfl"];

    // Derived data
    double dx = (xmax-xmin)/Nx;
    double dy = (ymax-ymin)/Ny;
    auto a = [](double x, double y) { return Eigen::Vector2d(y, -x); };
    double max_ax = ymax; // maximum of a_1 in domain
    double max_ay = -xmin; // maximum of a_2 in domain
    double dt_max = cfl / (max_ax/dx + max_ay/dy);
    std::string outfile = "u.txt";
    Eigen::MatrixXd u(Nx+2, Ny+2); // include ghost cells

    // apply initial condition to (1..Nx)x(1..Ny)
    //// ANCSE_START_TEMPLATE
    //----------------ICStart----------------
    for (int i = 1 ; i <= Nx ; i++) {
        for(int j = 1 ; j <= Ny ; j++) {
            u(i,j) = ic(xmin + (i-0.5)*dx, ymin + (j-0.5)*dy) ;
        }
    }
    //----------------ICEnd----------------
    //// ANCSE_END_TEMPLATE
    applyBoundaryConditions(u); // Complete u with BCs
    Eigen::MatrixXd u_old = u;

    wipeFile(outfile); // clear contents of output file

    double t = 0;
    std::vector<double> times;
    times.push_back(t);
    appendMatrixToFile(outfile, u.block(1,1,Nx,Ny));

    // Iterate over time
    while(t < T) {
        double dt = std::min(dt_max, T-t); // make sure we don't go beyond T
        t += dt;

        // Call updateUpwind. Don't forget the boundary conditions!
        //// ANCSE_START_TEMPLATE
        updateUpwind(u, u_old, dx, dy, dt, xmin, ymin, a);
        applyBoundaryConditions(u);
        //// ANCSE_END_TEMPLATE

        appendMatrixToFile(outfile, u.block(1,1,Nx,Ny));
        u_old = u;
        times.push_back(t);
    }
    writeToFile("time.txt", times);
}
