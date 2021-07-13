#include "Definitions.hpp"
#include "Matrix.hpp"
#include "Initialize.hpp"

namespace Initialize {

void initMovingwall(Matrix<Vec>& U, double uMax, size_t Nx, size_t Ny){
    /// Initialize Moving Wall
    for(size_t i = 0; i < Nx + 2; ++i){
         U(i, Ny + 1).comp[0] = uMax;
    }
}

// MICROSCOPIC INITIAL CONDITION
// TODO: use CS
void initDistfunc(Matrix<Node>& fIn, Matrix<Vec>& U, Matrix<double>& RHO,
                                                        size_t Nx, size_t Ny){
    double c_u;
    for (size_t i = 0; i < Nx + 2; ++i) {
        for (size_t j = 0; j < Ny + 2; ++j) {
            for (size_t q = 0; q < N_DIRECTIONS; ++q) {
                c_u = 3 * (U(i,j).comp[0] * LATTICE_VELOCITIES[q][0]
                        +  U(i,j).comp[1] * LATTICE_VELOCITIES[q][1]);

                fIn(i,j).dir[q] = RHO(i, j) * LATTICE_WEIGHTS[q]
                                   * ( 1 + c_u + 0.5 * (c_u * c_u)
                                      - 1.5 * ( U(i,j).comp[0] * U(i,j).comp[0] + U(i,j).comp[1] * U(i,j).comp[1] ));
            }
        }
    }
}

} // namespace Initialize
