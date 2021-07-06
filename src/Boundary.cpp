#include "Definitions.hpp"
#include "Matrix.hpp"
#include "Boundary.hpp"
#include "cuda.h"

namespace Boundary{
/// SETTING BOUNDARIES
// Moving wall
__global__ void setMovingwall(Node * fIn, Vec * U, double * RHO, double uMax, size_t Nx){

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    size_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if(j == 0) {
        // Macroscopic Dirichlet boundary conditions
        U[i].comp[0] = uMax;
        U[i].comp[1] = 0.0;
        RHO[i] = (fIn[i].dir[0] + fIn[i].dir[1] + fIn[i].dir[3]
        + 2*(fIn[i].dir[2] + fIn[i].dir[5] + fIn[i].dir[6])) / (1 - U[i].comp[1]);

        // Microscopic  Zou/He boundary conditions
        fIn[i].dir[4] = fIn[i].dir[2] - 2 / 3 * RHO[i] * U[i].comp[1];
        fIn[i].dir[7] = fIn[i].dir[5] + 0.5 * (fIn[i].dir[1] - fIn[i].dir[3])
                    - 0.5 * (RHO[i] * U[i].comp[0]) - 1/6 * (RHO[i] * U[i].comp[1]);
        fIn[i].dir[8] = fIn[i].dir[6] - 0.5 * (fIn[i].dir[1] - fIn[i].dir[3])
                    + 0.5 * (RHO[i] * U[i].comp[0]) - 1/6 * (RHO[i] * U[i].comp[1]);
    }
}


// Set bounce back cells
__global__ void setBounceback(Node * fOut, Node * fIn, size_t Nx, size_t Ny){

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    size_t j = blockIdx.y * blockDim.y + threadIdx.y;

    if (j == N + 1) {
        for(size_t q = 0; q < N_DIRECTIONS; ++q){
            fOut(i, Ny + 1).dir[q] = fIn(i, Ny + 1).dir[OPP[q]];
        }
    }

    if (i == 0) {
        for(size_t q = 0; q < N_DIRECTIONS; ++q){
            fOut(0, j).dir[q] = fIn(0, j).dir[OPP[q]];
        }
    }

    if (i == N + 1) {
        for(size_t q = 0; q < N_DIRECTIONS; ++q){
            fOut(Nx + 1, j).dir[q] = fIn(Nx + 1, j).dir[OPP[q]];
        }
    }
}

} // namespace Boundary
