#include "qclab/qgates/QGate1.hpp"
#include "apply.hpp"

namespace qclab::qgates {

  // apply
  template <typename T>
  void QGate1< T >::apply( Op op , const int nbQubits ,
                           std::vector< T >& vector , const int offset ) const {
    const int qubit = this->qubit() + offset ;
    auto f = lambda_QGate1( op , this->matrix() , vector.data() ) ;
    apply2( nbQubits , qubit , f ) ;
  }

#ifdef QCLAB_OMP_OFFLOADING
  // apply_device
  template <typename T>
  void QGate1< T >::apply_device( Op op , const int nbQubits , T* vector ,
                                  const int offset ) const {
    const int qubit = this->qubit() + offset ;
    auto f = lambda_QGate1( op , this->matrix() , vector ) ;
    apply_device2( nbQubits , qubit , f ) ;
  }
#endif

  // apply
  template <typename T>
  void QGate1< T >::apply( Side side , Op op , const int nbQubits ,
                           qclab::dense::SquareMatrix< T >& matrix ,
                           const int offset ) const {
    assert( nbQubits >= 1 ) ;
    assert( matrix.size() == 1 << nbQubits ) ;
    const int qubit = this->qubit() + offset ;
    assert( qubit < nbQubits ) ;
    // operation
    qclab::dense::SquareMatrix< T >  mat1 = this->matrix() ;
    qclab::dense::operateInPlace( op , mat1 ) ;
    // side
    if ( side == Side::Left ) {
      if ( nbQubits == 1 ) {
        matrix *= mat1 ;
      } else {
        // matrix *= kron( Ileft , mat1 , Iright )
        const int64_t nLeft  = 1 << qubit ;
        const int64_t nRight = 1 << ( nbQubits - qubit - 1 ) ;
        #pragma omp parallel for
        for ( int64_t i = 0; i < matrix.rows(); i++ ) {
          for ( int64_t k = 0; k < nLeft; k++ ) {
            for ( int64_t j = 0; j < nRight; j++ ) {
              const int64_t j1 = j + 2 * k * nRight ;
              const int64_t j2 = j1 + nRight ;
              const T x1 = matrix(i,j1) ;
              const T x2 = matrix(i,j2) ;
              matrix(i,j1) = mat1(0,0) * x1 + mat1(1,0) * x2 ;
              matrix(i,j2) = mat1(0,1) * x1 + mat1(1,1) * x2 ;
            }
          }
        }
      }
    } else {
      if ( nbQubits == 1 ) {
        matrix = mat1 * matrix ;
      } else {
        // matrix = kron( Ileft , mat1 , Iright ) * matrix
        const int64_t nLeft  = 1 << qubit ;
        const int64_t nRight = 1 << ( nbQubits - qubit - 1 ) ;
        #pragma omp parallel for
        for ( int64_t j = 0; j < matrix.cols(); j++ ) {
          for ( int64_t k = 0; k < nLeft; k++ ) {
            for ( int64_t i = 0; i < nRight; i++ ) {
              const int64_t i1 = i + 2 * k * nRight ;
              const int64_t i2 = i1 + nRight ;
              const T x1 = matrix(i1,j) ;
              const T x2 = matrix(i2,j) ;
              matrix(i1,j) = mat1(0,0) * x1 + mat1(0,1) * x2 ;
              matrix(i2,j) = mat1(1,0) * x1 + mat1(1,1) * x2 ;
            }
          }
        }
      }
    }
  }

  template class QGate1< float > ;
  template class QGate1< double > ;
  template class QGate1< std::complex< float > > ;
  template class QGate1< std::complex< double > > ;

} // namespace qclab::qgates

