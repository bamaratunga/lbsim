#ifndef __PROCESSES_HPP__
#define __PROCESSES_HPP__

#include "Definitions.hpp"
#include "Matrix.hpp"

namespace Processes{

    void computeDensity(const Node * currentNode, double * density);

    void computeVelocity(const Node * currentNode, Vec * velocity, double density);

    void computefEq(Node * eqField, const Vec * velocity, double density);

    void computePostCollisionDistributions(Node * outputNode, Node * currentNode, const Node * eqField, double omega);

    void calcQuantities(Matrix<Node>& fIn, Matrix<Vec>& U, Matrix<double>& RHO,
                                                            size_t Nx, size_t Ny);

    void doCollision(Matrix<Node>& fOut, Matrix<Node>& fIn, Matrix<Node>& fEq, Matrix<Vec>& U, Matrix<double>& RHO,
                                                                                  double omega, size_t Nx, size_t Ny);

    void doStreaming(Matrix<Node>& fIn, Matrix<Node>& fOut, size_t Nx, size_t Ny);

} // namespace Processes

#endif // __PROCESSES_HPP__
