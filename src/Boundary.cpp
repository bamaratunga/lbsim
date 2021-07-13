#include "Definitions.hpp"
#include "Matrix.hpp"
#include "Boundary.hpp"

namespace Boundary{

    /******************
     *   6   2   5
     *     \ | /
     *   3 - 0 - 1
     *     / | \
     *   7   4   8
    ******************/

/// SETTING BOUNDARIES
// Moving wall
void setMovingwall(Matrix<Node>& fIn, Matrix<Vec>& U, Matrix<double>& RHO, double uMax, size_t Nx, size_t Ny){
    for (size_t i = 0; i < Nx + 2; ++i) {
        // Macroscopic Dirichlet boundary conditions
        U(i, Ny + 1).comp[0] = uMax;
        U(i, Ny + 1).comp[1] = 0.0;
        RHO(i, Ny + 1) =   (fIn(i, Ny + 1).dir[0] + fIn(i, Ny + 1).dir[1] + fIn(i, Ny + 1).dir[3]
                  + 2*(fIn(i, Ny + 1).dir[2] + fIn(i, Ny + 1).dir[5] + fIn(i, Ny + 1).dir[6])) / (1 - U(i, Ny + 1).comp[1]);

        // Microscopic  Zou/He boundary conditions
        fIn(i, Ny + 1).dir[4] = fIn(i, Ny + 1).dir[2] - 2 / 3 * RHO(i, Ny + 1) * U(i, Ny + 1).comp[1];
        fIn(i, Ny + 1).dir[7] = fIn(i, Ny + 1).dir[5] + 0.5 * (fIn(i, Ny + 1).dir[1] - fIn(i, Ny + 1).dir[3])
                    - 0.5 * (RHO(i, Ny + 1) * U(i, Ny + 1).comp[0]) - 1/6 * (RHO(i, Ny + 1) * U(i, Ny + 1).comp[1]);
        fIn(i, Ny + 1).dir[8] = fIn(i, Ny + 1).dir[6] - 0.5 * (fIn(i, Ny + 1).dir[1] - fIn(i, Ny + 1).dir[3])
                    + 0.5 * (RHO(i, Ny + 1) * U(i, Ny + 1).comp[0]) - 1/6 * (RHO(i, Ny + 1) * U(i, Ny + 1).comp[1]);
    }
}


// Set bounce back cells
void setBounceback(Matrix<Node>& fOut, Matrix<Node>& fIn, size_t Nx, size_t Ny){
    for(size_t i = 0; i < Nx + 2; ++i){
        for(size_t q = 0; q < N_DIRECTIONS; ++q){
            fOut(i, 0).dir[q] = fIn(i, 0).dir[OPP[q]];
        }
    }
    for(size_t j = 0; j < Ny + 1; ++j){
        for(size_t q = 0; q < N_DIRECTIONS; ++q){
            fOut(0, j).dir[q] = fIn(0, j).dir[OPP[q]];
            fOut(Nx + 1, j).dir[q] = fIn(Nx + 1, j).dir[OPP[q]];
        }
    }
}

} // namespace Boundary
