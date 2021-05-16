//  (C) Copyright Roel Van Beeumen 2021.

#ifndef qclab_dense_kron_hpp
#define qclab_dense_kron_hpp

#include "SquareMatrix.hpp"
#include <cassert>

namespace qclab {

  namespace dense {

    /// Calculates the Kronecker product of `A` and `B`.
    template <typename T>
    void kron( const T& A , const T& B , T& kronAB ) {
      const int64_t rowsA = A.rows() ;
      const int64_t colsA = A.cols() ;
      const int64_t rowsB = B.rows() ;
      const int64_t colsB = B.cols() ;
      assert( kronAB.rows() == rowsA * rowsB ) ;
      assert( kronAB.cols() == colsA * colsB ) ;
      #pragma omp parallel for
      for ( int64_t i = 0; i < rowsA; i++ ) {        // i loops over rows of A
        for ( int64_t k = 0; k < rowsB; k++ ) {      // k loops over rows of B
          for ( int64_t j = 0; j < colsA; j++ ) {    // j loops over cols of A
            for ( int64_t l = 0; l < colsB; l++ ) {  // l loops over cols of B
              kronAB( i*rowsB + k , j*colsB + l ) = A(i,j) * B(k,l) ;
            }
          }
        }
      }
    } // kron(A,B,kronAB)

    /// Calculates the Kronecker product of `A` and `B`.
    template <typename T>
    T kron( const T& A , const T& B ) {
      T kronAB( A.rows() * B.rows() , A.cols() * B.cols() ) ;
      kron( A , B , kronAB ) ;
      return kronAB ;
    } // kron(A,B)

    /// Calculates the Kronecker product of `A` and `B`.
    template <typename T>
    SquareMatrix< T > kron( const SquareMatrix< T >& A ,
                            const SquareMatrix< T >& B ) {
      SquareMatrix< T > kronAB( A.rows() * B.rows() ) ;
      kron( A , B , kronAB ) ;
      return kronAB ;
    } // kron(A,B)

  } // namespace dense

} // namespace qclab

#endif

