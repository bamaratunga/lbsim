#include "Definitions.hpp"
#include "Matrix.hpp"
#include "Boundary.cpp"

namespace Boundary{
/// SETTING BOUNDARIES
// Moving wall
void setMovingwall(Matrix<Node> * fIn, Matrix<Vec> * U, Matrix<double> * RHO,
                                          double uMax, size_t Nx, size_t Ny){
    for (size_t i = 0; i < Nx + 2; ++i) {
        // Macroscopic Dirichlet boundary conditions
        U(i, 0).comp[0] = uMax;
        U(i, 0).comp[1] = 0.0;
        RHO(i, 0) =   (fIn(i, 0).dir[0] + fIn(i, 0).dir[1] + fIn(i, 0).dir[3]
                  + 2*(fIn(i, 0).dir[2] + fIn(i, 0).dir[5] + fIn(i, 0).dir[6])) / (1 - U(i, 0).comp[1]);

        // Microscopic  Zou/He boundary conditions
        fIn(i, 0).dir[4] = fIn(i, 0).dir[2] - 2 / 3 * RHO(i, 0) * U(i, 0).comp[1];
        fIn(i, 0).dir[7] = fIn(i, 0).dir[5] + 0.5 * (fIn(i, 0).dir[1] - fIn(i, 0).dir[3])
                    - 0.5 * (RHO(i, 0) * U(i, 0).comp[0]) - 1/6 * (RHO(i, 0) * U(i, 0).comp[1]);
        fIn(i, 0).dir[8] = fIn(i, 0).dir[6] + 0.5 * (fIn(i, 0).dir[1] - fIn(i, 0).dir[3])
                    + 0.5 * (RHO(i, 0) * U(i, 0).comp[0]) - 1/6 * (RHO(i, 0) * U(i, 0).comp[1]);
    }
}


// Set bounce back cells
void setBounceback(Matrix<Node> * fIn, Matrix<Node> * fOut, size_t Nx, size_t Ny){
    for(size_t i = 0; i < Nx + 2; ++i){
        for(size_t q = 0; q < N_DIRECTIONS; ++q){
            fOut(i, Ny + 1).dir[q] = fIn(i, Ny + 1).dir[OPP[q]];
        }
    }
    for(size_t j = 1; j < Ny + 2; ++j){
        for(size_t q = 0; q < N_DIRECTIONS; ++q){
            fOut(0, j).dir[q] = fIn(0, j).dir[OPP[q]];
            fOut(Nx + 1, j).dir[q] = fIn(Nx + 1, j).dir[OPP[q]];
        }
    }
}

} // namespace Boundary
