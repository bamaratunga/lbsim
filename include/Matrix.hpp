#ifndef __MATRIX_HPP__
#define __MATRIX_HPP__

#include <cstdlib>

template <typename T> class Matrix {

  public:
    Matrix<T>() = default;

    /// Constructor without an initial value.
    Matrix<T>(size_t i_max, size_t j_max) : _imax(i_max), _jmax(j_max) {
        _container = (T *)malloc(_imax * _jmax * sizeof(T));
    }

    /// Element access and modify using index
    T &operator()(size_t i, size_t j) {
        return _container[j * _imax + i];
    }

    /// Access of the size of the structure
    size_t size() const { return _imax * _jmax; }

    /// get the number of elements in x direction
    size_t imax() const { return _imax; }

    /// get the number of elements in y direction
    size_t jmax() const { return _jmax; }

  private:
    /// Number of elements in x direction
    size_t _imax;
    /// Number of elements in y direction
    size_t _jmax;

    /// Data container
    T * _container;
};

#endif // __MATRIX_HPP__
