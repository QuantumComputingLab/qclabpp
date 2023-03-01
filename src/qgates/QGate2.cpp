#include "qclab/qgates/QGate2.hpp"
#include "apply.hpp"

namespace qclab::qgates {

  // apply
  template <typename T>
  void QGate2< T >::apply( Op op , const int nbQubits ,
                           std::vector< T >& vector , const int offset ) const {
    auto qubits = this->qubits() ;
    qubits[0] += offset ;
    qubits[1] += offset ;
    auto f = lambda_QGate2( op , this->matrix() , vector.data() ) ;
    apply4( nbQubits , qubits[0] , qubits[1] , f ) ;
  }

#ifdef QCLAB_OMP_OFFLOADING
  // apply_device
  template <typename T>
  void QGate2< T >::apply_device( Op op , const int nbQubits , T* vector ,
                                  const int offset ) const {
    auto qubits = this->qubits() ;
    qubits[0] += offset ;
    qubits[1] += offset ;
    auto f = lambda_QGate2( op , this->matrix() , vector ) ;
    apply_device4( nbQubits , qubits[0] , qubits[1] , f ) ;
  }
#endif

  // apply
  template <typename T>
  void QGate2< T >::apply( Side side , Op op , const int nbQubits ,
                           qclab::dense::SquareMatrix< T >& matrix ,
                           const int offset ) const {
    assert( nbQubits >= 2 ) ;
    assert( matrix.size() == 1 << nbQubits ) ;
    auto qubits = this->qubits() ;
    qubits[0] += offset ;
    qubits[1] += offset ;
    assert( qubits[0] < nbQubits ) ; assert( qubits[1] < nbQubits ) ;
    assert( qubits[0] + 1 == qubits[1] ) ;  // nearest neighbor qubtis
    // operation
    qclab::dense::SquareMatrix< T >  mat2 = this->matrix() ;
    qclab::dense::operateInPlace( op , mat2 ) ;
    // side
    if ( side == Side::Left ) {
      if ( nbQubits == 2 ) {
        matrix *= mat2 ;
      } else {
        // matrix *= kron( Ileft , mat2 , Iright )
        const int64_t nLeft  = 1 << qubits[0] ;
        const int64_t nRight = 1 << ( nbQubits - qubits[1] - 1 ) ;
        #pragma omp parallel for
        for ( int64_t i = 0; i < matrix.rows(); i++ ) {
          for ( int64_t k = 0; k < nLeft; k++ ) {
            for ( int64_t j = 0; j < nRight; j++ ) {
              const int64_t j1 = j + 4 * k * nRight ;
              const int64_t j2 = j1 + nRight ;
              const int64_t j3 = j2 + nRight ;
              const int64_t j4 = j3 + nRight ;
              const T x1 = matrix(i,j1) ;
              const T x2 = matrix(i,j2) ;
              const T x3 = matrix(i,j3) ;
              const T x4 = matrix(i,j4) ;
              matrix(i,j1) = mat2(0,0) * x1 + mat2(1,0) * x2 +
                             mat2(2,0) * x3 + mat2(3,0) * x4 ;
              matrix(i,j2) = mat2(0,1) * x1 + mat2(1,1) * x2 +
                             mat2(2,1) * x3 + mat2(3,1) * x4 ;
              matrix(i,j3) = mat2(0,2) * x1 + mat2(1,2) * x2 +
                             mat2(2,2) * x3 + mat2(3,2) * x4 ;
              matrix(i,j4) = mat2(0,3) * x1 + mat2(1,3) * x2 +
                             mat2(2,3) * x3 + mat2(3,3) * x4 ;
            }
          }
        }
      }
    } else {
      if ( nbQubits == 2 ) {
        matrix = mat2 * matrix ;
      } else {
        // matrix = kron( Ileft , mat2 , Iright ) * matrix
        const int64_t nLeft  = 1 << qubits[0] ;
        const int64_t nRight = 1 << ( nbQubits - qubits[1] - 1 ) ;
        #pragma omp parallel for
        for ( int64_t j = 0; j < matrix.cols(); j++ ) {
          for ( int64_t k = 0; k < nLeft; k++ ) {
            for ( int64_t i = 0; i < nRight; i++ ) {
              const int64_t i1 = i + 4 * k * nRight ;
              const int64_t i2 = i1 + nRight ;
              const int64_t i3 = i2 + nRight ;
              const int64_t i4 = i3 + nRight ;
              const T x1 = matrix(i1,j) ;
              const T x2 = matrix(i2,j) ;
              const T x3 = matrix(i3,j) ;
              const T x4 = matrix(i4,j) ;
              matrix(i1,j) = mat2(0,0) * x1 + mat2(0,1) * x2 +
                             mat2(0,2) * x3 + mat2(0,3) * x4 ;
              matrix(i2,j) = mat2(1,0) * x1 + mat2(1,1) * x2 +
                             mat2(1,2) * x3 + mat2(1,3) * x4 ;
              matrix(i3,j) = mat2(2,0) * x1 + mat2(2,1) * x2 +
                             mat2(2,2) * x3 + mat2(2,3) * x4 ;
              matrix(i4,j) = mat2(3,0) * x1 + mat2(3,1) * x2 +
                             mat2(3,2) * x3 + mat2(3,3) * x4 ;
            }
          }
        }
      }
    }
  }

  template class QGate2< float > ;
  template class QGate2< double > ;
  template class QGate2< std::complex< float > > ;
  template class QGate2< std::complex< double > > ;

} // namespace qclab::qgates

