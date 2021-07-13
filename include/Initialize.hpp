#ifndef __INITIALIZE_HPP__
#define __INITIALIZE_HPP__

#include "Definitions.hpp"
#include "Matrix.hpp"

namespace Initialize {

    void initMovingwall(Matrix<Vec>& U, double uMax, size_t Nx, size_t Ny);

    void initDistfunc(Matrix<Node>& fIn, Matrix<Vec>& U, Matrix<double>& RHO,
                                                            size_t Nx, size_t Ny);

} // namespace Initialzie

#endif // __INITIALIZE_HPP__
