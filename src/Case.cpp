#include <iostream>
#include <cstdlib>
#include <cmath>

#include "Definitions.hpp"
#include "Matrix.hpp"
#include "Utils.hpp"
#include "Initialize.hpp"
#include "Boundary.hpp"
#include "Processes.hpp"


int main(){

    /// DIMENSIONS
    size_t Nx      = 50;         // number of cells in x-direction
    size_t Ny      = 50;         // number of cells in y-direction
    size_t Nt      = 2000;        // number of time steps
    size_t plot_interval = 100;   // How many timesteps for the next plot update

    /// FLOW PARAMETERS
    double Re      = 300;               // Reynolds number (Main flow parameter)0000
    double uMax    = 0.08;              // Velocity of lid
    double nu      = uMax * Nx / Re;    // kinematic viscosity (is chosen such that the given Raynolds number is achieved)
    double omega   = 2 / (6 * nu + 1);  // relaxation parameter

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
    Initialize::initMovingwall(U, uMax, Nx);
    // Microscopic quantities
    Initialize::initDistfunc(fIn, U, RHO, Nx, Ny);

/*************************************************************************************/

/*************************************************************************************/

    for(size_t t_step = 0; t_step < Nt; ++t_step){

        /// CALCULATE MACROSCOPIC QUANTITIES FOR EACH CELL
        Processes::calcQuantities(fIn, U, RHO, Nx, Ny);

        /// SETTING BOUNDARIES
        // Moving wall
        Boundary::setMovingwall(fIn, U, RHO, uMax, Nx);

        /// COLLISION STEP
        Processes::doCollision(fEq, fOut, U, RHO, omega, Nx, Ny);

        /// SETTING BOUNDARIES
        // Set bounce back cells
        Boundary::setBounceback(fIn, fOut, Nx, Ny);

        // STREAMING STEP
        Processes::doStreaming(fIn, fOut, Nx, Ny);

        if(t_step % plot_interval == 0){
            Utils::writeVtkOutput(U, Nx, Ny, t_step);
        }
    }
}
