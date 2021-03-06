#include <iostream>
#include <cstdlib>
#include <chrono>
#include <cmath>

#include "Definitions.hpp"
#include "Matrix.hpp"
#include "Utils.hpp"
#include "Initialize.hpp"
#include "Boundary.hpp"
#include "Processes.hpp"


int main(int argn, char **args){

    if (argn <= 1) {
        std::cout << "ERROR: No input file is provided. Exiting..." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string file_name{args[1]};
    Inputs input;
    Utils::read_inputs(file_name, input);

    /// DIMENSIONS
    size_t Nx      = input.x_dim;        // number of cells in x-direction
    size_t Ny      = input.y_dim;        // number of cells in y-direction
    size_t Nt      = input.timesteps;    // number of time steps

    bool output_results = input.output_results;
    size_t plot_interval = input.plot_interval;   // How many timesteps for the next plot update

    /// FLOW PARAMETERS
    double Re      = input.reynolds;    // Reynolds number (Main flow parameter)
    double uMax    = input.wall_vel;    // Velocity of lid
    double nu      = uMax * Nx / Re;    // kinematic viscosity (is chosen such that the given Raynolds number is achieved)
    double omega   = 2 / (6 * nu + 1);  // relaxation parameter

    auto start = std::chrono::high_resolution_clock::now();

    /// Input distribution:
    Matrix<Node> fIn  = Matrix<Node>(Nx + 2, Ny + 2);
    /// Output distributions (post collision):
    Matrix<Node> fOut = Matrix<Node>(Nx + 2, Ny + 2);
    /// Equilibrium distribution:
    Matrix<Node> fEq  = Matrix<Node>(Nx + 2, Ny + 2);

    /// Macroscopic velocity vector
    Matrix<Vec> U = Matrix<Vec>(Nx + 2, Ny + 2);
    /// Density
    Matrix<double> RHO = Matrix<double>(Nx + 2, Ny + 2);

    /// Zero initialize all
    for (size_t i = 0; i < Nx + 2; ++i) {
        for (size_t j = 0; j < Ny + 2; ++j) {
            for(size_t q = 0; q < N_DIRECTIONS; ++q){
                fIn(i, j).dir[q] = 0.0;
                fOut(i, j).dir[q] = 0.0;
                fEq(i, j).dir[q] = 0.0;
            }
            for(size_t d = 0; d < N_DIM; ++d){
                U(i, j).comp[d] = 0.0;
            }
            RHO(i, j) = 1.0;
        }
    }

    /// SET INITIAL CONDITIONS
    // Moving wall velocity
    Initialize::initMovingwall(U, uMax, Nx, Ny);
    // Microscopic quantities
    Initialize::initDistfunc(fIn, U, RHO, Nx, Ny);

/*************************************************************************************
*   SIMULATION
*************************************************************************************/

    std::cout << "Running simulation..." << std::endl;
    for(size_t t_step = 0; t_step < Nt; ++t_step){

        /// CALCULATE MACROSCOPIC QUANTITIES FOR EACH CELL
        Processes::calcQuantities(fIn, U, RHO, Nx, Ny);

        /// SETTING BOUNDARIES
        // Moving wall
        Boundary::setMovingwall(fIn, U, RHO, uMax, Nx, Ny);

        /// COLLISION STEP
        Processes::doCollision(fOut, fIn, fEq, U, RHO, omega, Nx, Ny);

        /// SETTING BOUNDARIES
        // Set bounce back cells
        Boundary::setBounceback(fOut, fIn, Nx, Ny);

        // STREAMING STEP
        Processes::doStreaming(fIn, fOut, Nx, Ny);

        if(output_results && t_step % plot_interval == 0){
            std::cout << "Writing vtk file at t = " << t_step << std::endl;
            Utils::writeVtkOutput(U, Nx, Ny, t_step, 0, input.case_name, input.dict_name);
        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << std::endl;
    std::cout << "LBSIM CPP VERSION" << std::endl;
    std::cout << "Simulating " << input.case_name << std::endl;
    std::cout << "Dimensions = " << Nx + 2 << " x " << Ny + 2 << std::endl;
    std::cout << "No. of iterations = " << Nt << std::endl;
    std::cout << "Execution time = " << duration.count() / 1e6 << "s" << std::endl;
}
