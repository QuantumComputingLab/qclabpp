//  (C) Copyright Roel Van Beeumen 2021.

#ifndef qclab_dense_transpose_hpp
#define qclab_dense_transpose_hpp

#include "SquareMatrix.hpp"
#include "util.hpp"
#include <complex>

namespace qclab {

  namespace dense {

    /// Performs the operation `op` in-place on the matrix `A`.
    template <typename T>
    void operateInPlace( Op op , T& A ) {
      if ( op == Op::NoTrans ) return ;
      const int64_t rows = A.rows() ;
      const int64_t cols = A.cols() ;
      if constexpr ( qclab::is_complex_v< typename T::value_type > ) {
        if ( op == Op::ConjTrans ) {
          // conjugate
          #pragma omp parallel for
          for ( int64_t j = 0; j < cols; j++ ) {
            for ( int64_t i = 0; i < rows; i++ ) {
              A(i,j) = std::conj( A(i,j) ) ;
            }
          }
        }
      }
      // transpose
      #pragma omp parallel for
      for ( int64_t j = 0; j < cols - 1; j++ ) {
        for ( int64_t i = j + 1; i < rows; i++ ) {
          std::swap( A(i,j) , A(j,i) ) ;
        }
      }
    }

    /// Performs the operation `op` on the matrix `A`.
    template <typename T>
    T operate( Op op , const T& A ) {
      T Aop( A ) ;
      operateInPlace( op , Aop ) ;
      return Aop ;
    }

    /// Calculates in-place the transpose of the matrix `A`.
    template <typename T>
    void transInPlace( T& A ) {
      operateInPlace( Op::Trans , A ) ;
    }

    /// Calculates transpose of the matrix `A`.
    template <typename T>
    T trans( const T& A ) {
      T Atrans( A ) ;
      operateInPlace( Op::Trans , Atrans ) ;
      return Atrans ;
    }

    /// Calculates in-place the conjugate transpose of the matrix`A`.
    template <typename T>
    void conjTransInPlace( T& A ) {
      operateInPlace( Op::ConjTrans , A ) ;
    }

    /// Calculates the conjugate transpose of the matrix`A`.
    template <typename T>
    T conjTrans( const T& A ) {
      T AconjTrans( A ) ;
      operateInPlace( Op::ConjTrans , AconjTrans ) ;
      return AconjTrans ;
    }

  } // namespace dense

} // namespace qclab

#endif

