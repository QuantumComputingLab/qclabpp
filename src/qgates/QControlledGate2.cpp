#include "qclab/qgates/QControlledGate2.hpp"
#include "apply.hpp"

namespace qclab::qgates {

  // apply
  template <typename T>
  void QControlledGate2< T >::apply( Op op ,
                                     const int nbQubits ,
                                     std::vector< T >& vector ,
                                     const int offset ) const {
    auto qubits = this->qubits() ;
    qubits[0] += offset ;
    qubits[1] += offset ;
    const int control = this->control() + offset ;
    const int target  = this->target()  + offset ;
    auto f = lambda_QGate1( op , this->gate()->matrix() , vector.data() ) ;
    apply4( nbQubits , qubits[0] , qubits[1] , control , target ,
            this->controlState() , f ) ;
  }

#ifdef QCLAB_OMP_OFFLOADING
  // apply_device
  template <typename T>
  void QControlledGate2< T >::apply_device( Op op ,
                                            const int nbQubits ,
                                            T* vector ,
                                            const int offset ) const {
    auto qubits = this->qubits() ;
    qubits[0] += offset ;
    qubits[1] += offset ;
    const int control = this->control() + offset ;
    const int target  = this->target()  + offset ;
    auto f = lambda_QGate1( op , this->gate()->matrix() , vector ) ;
    apply_device4( nbQubits , qubits[0] , qubits[1] , control , target ,
                   this->controlState() , f ) ;
  }
#endif

  // apply
  template <typename T>
  void QControlledGate2< T >::apply( Side side , Op op , const int nbQubits ,
                                     qclab::dense::SquareMatrix< T >& matrix ,
                                     const int offset ) const {
    using matrix_type = qclab::dense::SquareMatrix< T > ;
    assert( nbQubits >= 2 ) ;
    assert( matrix.size() == 1 << nbQubits ) ;
    auto qubits = this->qubits() ;
    qubits[0] += offset ;
    qubits[1] += offset ;
    assert( qubits[0] < nbQubits ) ; assert( qubits[1] < nbQubits ) ;
    // operation
    matrix_type  mat1 = this->gate()->matrix() ;
    qclab::dense::operateInPlace( op , mat1 ) ;
    // control / target
    if ( control() < target() ) {
      // control() < target()
      const int64_t d = 1 << ( target() - control() ) ;
      const int64_t nLeft  = 1 << target() ;
      const int64_t nRight = 1 << ( nbQubits - target() - 1 ) ;
      if ( side == Side::Left ) {
        #pragma omp parallel for
        for ( int64_t i = 0; i < matrix.rows(); i++ ) {
          for ( int64_t k = 0; k < nLeft; k++ ) {
            if ( controlState_ == 0 && ( k % d >= d/2 ) ) continue ;
            if ( controlState_ == 1 && ( k % d <  d/2 ) ) continue ;
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
      } else {
        #pragma omp parallel for
        for ( int64_t j = 0; j < matrix.cols(); j++ ) {
          for ( int64_t k = 0; k < nLeft; k++ ) {
            if ( controlState_ == 0 && ( k % d >= d/2 ) ) continue ;
            if ( controlState_ == 1 && ( k % d <  d/2 ) ) continue ;
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
    } else {
      // control() >= target()
      const int64_t d = 1 << ( nbQubits - control() ) ;
      const int64_t nLeft  = 1 << target() ;
      const int64_t nRight = 1 << ( nbQubits - target() - 1 ) ;
      if ( side == Side::Left ) {
        #pragma omp parallel for
        for ( int64_t i = 0; i < matrix.cols(); i++ ) {
          for ( int64_t k = 0; k < nLeft; k++ ) {
            for ( int64_t j = 0; j < nRight; j++ ) {
              if ( controlState_ == 0 && ( j % d >= d/2 ) ) continue ;
              if ( controlState_ == 1 && ( j % d <  d/2 ) ) continue ;
              const int64_t j1 = j + 2 * k * nRight ;
              const int64_t j2 = j1 + nRight ;
              const T x1 = matrix(i,j1) ;
              const T x2 = matrix(i,j2) ;
              matrix(i,j1) = mat1(0,0) * x1 + mat1(1,0) * x2 ;
              matrix(i,j2) = mat1(0,1) * x1 + mat1(1,1) * x2 ;
            }
          }
        }
      } else {
        #pragma omp parallel for
        for ( int64_t j = 0; j < matrix.cols(); j++ ) {
          for ( int64_t k = 0; k < nLeft; k++ ) {
            for ( int64_t i = 0; i < nRight; i++ ) {
              if ( controlState_ == 0 && ( i % d >= d/2 ) ) continue ;
              if ( controlState_ == 1 && ( i % d <  d/2 ) ) continue ;
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

  template class QControlledGate2< float > ;
  template class QControlledGate2< double > ;
  template class QControlledGate2< std::complex< float > > ;
  template class QControlledGate2< std::complex< double > > ;

} // namespace qclab::qgates

