#include "qclab/qgates/Phase45.hpp"
#include "apply.hpp"

namespace qclab::qgates {

  // apply
  template <typename T>
  void Phase45< T >::apply( Op op , const int nbQubits ,
                            std::vector< T >& vector ,
                            const int offset ) const {
    const int qubit = this->qubit() + offset ;
    auto f = lambda_Phase45( op , vector.data() ) ;
    apply2b( nbQubits , qubit , f ) ;
  }

#ifdef QCLAB_OMP_OFFLOADING
  // apply_device
  template <typename T>
  void Phase45< T >::apply_device( Op op , const int nbQubits , T* vector ,
                                   const int offset ) const {
    const int qubit = this->qubit() + offset ;
    auto f = lambda_Phase45( op , vector ) ;
    apply_device2b( nbQubits , qubit , f ) ;
  }
#endif

  // apply
  template <typename T>
  void Phase45< T >::apply( Side side , Op op , const int nbQubits ,
                            qclab::dense::SquareMatrix< T >& matrix ,
                            const int offset ) const {
    QGate1< T >::apply( side , op , nbQubits , matrix , offset ) ;
  }

  template class Phase45< std::complex< float > > ;
  template class Phase45< std::complex< double > > ;

} // namespace qclab::qgates

