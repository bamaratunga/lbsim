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

__device__ void computefEq(Node * eqField, const Vec * velocity, double density){

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

__device__ void computePostCollisionDistributions(Node * outputNode, Node * currentNode, const Node * eqField, double omega){

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

__global__ void doCollision(Node * fOut, Node * fIn, Node * fEq, Vec * U, double * RHO,
                                                            double omega, size_t Nx, size_t Ny) {
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    size_t j = blockIdx.y * blockDim.y + threadIdx.y;
    /// COLLISION STEP
    computefEq( fEq[i + j * (Nx + 2)], U[i + j * (Nx + 2)], RHO[i + j * (Nx + 2)] );
    computePostCollisionDistributions( fOut[i + j * (Nx + 2)], fIn[i + j * (Nx + 2)], fEq[i + j * (Nx + 2)], omega );
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
