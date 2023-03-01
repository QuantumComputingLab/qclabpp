#include "qclab/qgates/Phase.hpp"
#include "apply.hpp"

namespace qclab::qgates {

  // apply
  template <typename T>
  void Phase< T >::apply( Op op , const int nbQubits ,
                          std::vector< T >& vector , const int offset ) const {
    const int qubit = this->qubit() + offset ;
    T lambda = T( cos() , sin() ) ;
    auto f = lambda_Phase( op , lambda , vector.data() ) ;
    apply2b( nbQubits , qubit , f ) ;
  }

#ifdef QCLAB_OMP_OFFLOADING
  // apply_device
  template <typename T>
  void Phase< T >::apply_device( Op op , const int nbQubits , T* vector ,
                                 const int offset ) const {
    const int qubit = this->qubit() + offset ;
    T lambda = T( cos() , sin() ) ;
    auto f = lambda_Phase( op , lambda , vector ) ;
    apply_device2b( nbQubits , qubit , f ) ;
  }
#endif

  // apply
  template <typename T>
  void Phase< T >::apply( Side side , Op op , const int nbQubits ,
                          qclab::dense::SquareMatrix< T >& matrix ,
                          const int offset ) const {
    QGate1< T >::apply( side , op , nbQubits , matrix , offset ) ;
  }

  template class Phase< std::complex< float > > ;
  template class Phase< std::complex< double > > ;

} // namespace qclab::qgates

