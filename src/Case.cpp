#include <iostream>
#include <cstdlib>
#include <cmath>

#include "Definitions.hpp"
#include "Matrix.hpp"
#include "Utils.hpp"
#include "Initialize.hpp"
#include "Boundary.hpp"
#include "Processes.hpp"
#include "cuda.h"

#define gpuErrChk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess)
  {
    fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
    if (abort) exit(code);
  }
}

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
    size_t plot_interval = input.plot_interval;   // How many timesteps for the next plot update

    /// FLOW PARAMETERS
    double Re      = input.reynolds;    // Reynolds number (Main flow parameter)
    double uMax    = input.wall_vel;    // Velocity of lid
    double nu      = uMax * Nx / Re;    // kinematic viscosity (is chosen such that the given Raynolds number is achieved)
    double omega   = 2 / (6 * nu + 1);  // relaxation parameter

    /// Input distribution:
    Node * fIn;
    gpuErrChk(cudaMallocManaged(&fIn, (Nx + 2) * (Ny + 2) * sizeof(Node)));
    /// Output distributions (post collision):
    Node * fOut;
    gpuErrChk(cudaMallocManaged(&fOut, (Nx + 2) * (Ny + 2) * sizeof(Node)));
    /// Equilibrium distribution:
    Node * fEq;
    gpuErrChk(cudaMallocManaged(&fEq, (Nx + 2) * (Ny + 2) * sizeof(Node)));

    /// Macroscopic velocity vector
    Vec * U;
    gpuErrChk(cudaMallocManaged(&U, (Nx + 2) * (Ny + 2) * sizeof(Vec)));
    /// Density
    double * RHO;
    gpuErrChk(cudaMallocManaged(&RHO, (Nx + 2) * (Ny + 2) * sizeof(double)));

    /// Zero initialize all
    for (size_t i = 0; i < Nx + 2; ++i) {
        for (size_t j = 0; j < Ny + 2; ++j) {
            for(size_t q = 0; q < N_DIRECTIONS; ++q){
                fIn[i + j * (Nx + 2)].dir[q] = 0.0;
                fOut[i + j * (Nx + 2)].dir[q] = 0.0;
                fEq[i + j * (Nx + 2)].dir[q] = 0.0;
            }
            for(size_t d = 0; d < N_DIM; ++d){
                U[i + j * (Nx + 2)].comp[d] = 0.0;
            }
            RHO[i + j * (Nx + 2)] = 1.0;
        }
    }

    /// SET INITIAL CONDITIONS
    // Moving wall velocity
    Initialize::initMovingwall(U, uMax, Nx);
    // Microscopic quantities
    Initialize::initDistfunc(fIn, U, RHO, Nx, Ny);

/*************************************************************************************
*   SIMULATION
*************************************************************************************/

    for(size_t t_step = 0; t_step < Nt; ++t_step){

        /// CALCULATE MACROSCOPIC QUANTITIES FOR EACH CELL
        Processes::calcQuantities(fIn, U, RHO, Nx, Ny);

        /// SETTING BOUNDARIES
        // Moving wall
        Boundary::setMovingwall(fIn, U, RHO, uMax, Nx);

        /// COLLISION STEP
        Processes::doCollision(fOut, fIn, fEq, U, RHO, omega, Nx, Ny);

        /// SETTING BOUNDARIES
        // Set bounce back cells
        Boundary::setBounceback(fOut, fIn, Nx, Ny);

        // STREAMING STEP
        Processes::doStreaming(fIn, fOut, Nx, Ny);

        if(t_step % plot_interval == 0){
            std::cout << "Writing vtk file at t = " << t_step << std::endl;
            Utils::writeVtkOutput(U, Nx, Ny, t_step, 0, input.case_name, input.dict_name);
        }
    }
}
