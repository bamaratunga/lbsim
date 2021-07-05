#include "Definitions.hpp"
#include "Matrix.hpp"
#include "Processes.hpp"

namespace Processes {

void computeDensity(const Node * currentNode, double * density){
    *density = 0;
    for (int q = 0; q < N_DIRECTIONS; ++q) {
        *density += currentNode->dir[q];
    }
}

void computeVelocity(const Node * currentNode, Vec * velocity, double density) {

    for (int d = 0; d < N_DIM; ++d) {
        // Initialize each velocity to zero
        velocity->comp[d] = 0.0;
        for (int q = 0; q < N_DIRECTIONS; ++q) {
           velocity->comp[d] += currentNode->dir[q] * LATTICE_VELOCITIES[q][d];
        }

        velocity->comp[d] = (density == 0) ? 0: velocity->comp[d] / density;
    }
}

void computefEq(Node * eqField, const Vec * velocity, double density){

	for (int q = 0; q < N_DIRECTIONS; ++q) {

        double c_dot_u = 0.0;
        double u_dot_u = 0.0;
		for (int d = 0; d < N_DIM; ++d) {
            c_dot_u += LATTICE_VELOCITIES[q][d] * velocity->comp[d];
            u_dot_u += velocity->comp[d] * velocity->comp[d];
		}

		eqField->dir[q] =  LATTICE_WEIGHTS[q] * density * (1 + c_dot_u / (CS * CS)
                    + (c_dot_u * c_dot_u)/(2 * CS * CS * CS * CS)
                    - u_dot_u / (2 * CS * CS));
	}
}

void computePostCollisionDistributions(Node * outputNode, Node * currentNode, const Node * eqField, double omega){

    for (int q = 0; q < N_DIRECTIONS; ++q) {
        outputNode->dir[q] = currentNode->dir[q] - omega *( currentNode->dir[q] - eqField->dir[q] );
    }
}

void calcQuantities(Matrix<Node>& fIn, Matrix<Vec>& U, Matrix<double>& RHO,
                                                        size_t Nx, size_t Ny){
    /// CALCULATE MACROSCOPIC QUANTITIES FOR EACH CELL
    for (size_t i = 0; i < Nx + 2; ++i) {
        for (size_t j = 0; j < Ny + 2; ++j) {
            computeDensity( &fIn(i, j), &RHO(i, j));
            computeVelocity( &fIn(i, j), &U(i, j), RHO(i, j));
        }
    }
}

void doCollision(Matrix<Node>& fOut, Matrix<Node>& fIn, Matrix<Node>& fEq, Matrix<Vec>& U, Matrix<double>& RHO,
                                                           double omega, size_t Nx, size_t Ny) {
    /// COLLISION STEP
    for (size_t i = 0; i < Nx + 2; ++i) {
        for (size_t j = 0; j < Ny + 2; ++j) {
            computefEq( &fEq(i, j), &U(i, j), RHO(i, j) );
            computePostCollisionDistributions( &fOut(i, j), &fIn(i, j), &fEq(i,j), omega );
        }
    }
}

void doStreaming(Matrix<Node>& fOut, Matrix<Node>& fIn, size_t Nx, size_t Ny){

    int Cx, Cy;
    for (size_t j = 0; j < Ny + 2; ++j) {
        for (size_t i = 0; i < Nx + 2; ++i) {
            for (size_t q = 0; q < N_DIRECTIONS; ++q) {
                Cx = LATTICE_VELOCITIES[q][0];
                Cy = LATTICE_VELOCITIES[q][1];
                fOut(i, j).dir[q] = fIn( (i - Cx + Nx + 2) % (Nx + 2), (j - Cy + Ny + 2) % (Ny + 2) ).dir[q];
            }
        }
    }
}

} // namespace Processes
