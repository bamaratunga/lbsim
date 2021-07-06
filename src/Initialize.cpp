#include "Definitions.hpp"
#include "Matrix.hpp"
#include "Initialize.hpp"
#include "cuda.h"

namespace Initialize {

__host__ void initMovingwall(Vec * U, double uMax, size_t Nx){
    /// Initialize Moving Wall
    for(size_t i = 0; i < Nx + 2; ++i){
         U[i + 0 * (Nx + 2)].comp[0] = uMax;
    }
}

// MICROSCOPIC INITIAL CONDITION
// TODO: use CS
__global__ void initDistfunc(Node * fIn, Vec * U, double RHO, size_t Nx, size_t Ny){

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    size_t j = blockIdx.y * blockDim.y + threadIdx.y;

    double c_u;
    for (size_t q = 0; q < N_DIRECTIONS; ++q) {
        c_u = (U[i + j * (Nx + 2)].comp[0] * LATTICE_VELOCITIES[q][0]
            +  U[i + j * (Nx + 2)].comp[1] * LATTICE_VELOCITIES[q][1]);

        fIn[i + j * (Nx + 2)].dir[q] = RHO[i + j * (Nx + 2)] * LATTICE_WEIGHTS[q]
                              * ( 1 + c_u / (CS * CS) + (c_u * c_u) / (2 * CS * CS * CS * CS)
                              - ( U[i + j * (Nx + 2)].comp[0] * U[i + j * (Nx + 2)].comp[0]
                              +   U[i + j * (Nx + 2)].comp[1] * U[i + j * (Nx + 2)].comp[1] ) / (2 * CS * CS));
    }
}

} // namespace Initialize
