#ifndef __BOUNDARY_HPP__
#define __BOUNDARY_HPP__

namespace Boundary{

    void setMovingwall(Matrix<Node> * fIn, Matrix<Vec> * U, Matrix<double> * RHO,
                                              double uMax, size_t Nx, size_t Ny);

    void setBounceback(Matrix<Node> * fIn, Matrix<Node> * fOut, size_t Nx, size_t Ny);

} // namespace Boundary

#endif // __BOUNDARY_HPP__
