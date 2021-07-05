#ifndef __BOUNDARY_HPP__
#define __BOUNDARY_HPP__

namespace Boundary{

    void setMovingwall(Matrix<Node>& fIn, Matrix<Vec>& U, Matrix<double>& RHO, double uMax, size_t Nx);

    void setBounceback(Matrix<Node>& fOut, Matrix<Node>& fIn, size_t Nx, size_t Ny);

} // namespace Boundary

#endif // __BOUNDARY_HPP__
