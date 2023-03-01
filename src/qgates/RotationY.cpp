#include "qclab/qgates/RotationY.hpp"
#include "apply.hpp"

namespace qclab::qgates {

  // apply
  template <typename T>
  void RotationY< T >::apply( Op op , const int nbQubits ,
                              std::vector< T >& vector ,
                              const int offset ) const {
    const int qubit = this->qubit() + offset ;
    auto f = lambda_RotationY( op , this->cos() , this->sin() , vector.data() );
    apply2( nbQubits , qubit , f ) ;
  }

#ifdef QCLAB_OMP_OFFLOADING
  // apply_device
  template <typename T>
  void RotationY< T >::apply_device( Op op , const int nbQubits , T* vector ,
                                     const int offset ) const {
    const int qubit = this->qubit() + offset ;
    auto f = lambda_RotationY( op , this->cos() , this->sin() , vector ) ;
    apply_device2( nbQubits , qubit , f ) ;
  }
#endif

  // apply
  template <typename T>
  void RotationY< T >::apply( Side side , Op op , const int nbQubits ,
                              qclab::dense::SquareMatrix< T >& matrix ,
                              const int offset ) const {
    QGate1< T >::apply( side , op , nbQubits , matrix , offset ) ;
  }

  template class RotationY< float > ;
  template class RotationY< double > ;
  template class RotationY< std::complex< float > > ;
  template class RotationY< std::complex< double > > ;

} // namespace qclab::qgates

