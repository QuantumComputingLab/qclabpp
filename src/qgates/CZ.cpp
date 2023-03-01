#include "qclab/qgates/CZ.hpp"
#include "apply.hpp"

namespace qclab::qgates {

  // apply
  template <typename T>
  void CZ< T >::apply( Op op , const int nbQubits , std::vector< T >& vector ,
                       const int offset ) const {
    auto qubits = this->qubits() ;
    qubits[0] += offset ;
    qubits[1] += offset ;
    const int control = this->control() + offset ;
    const int target  = this->target()  + offset ;
    auto f = lambda_PauliZ( op , vector.data() ) ;
    apply4b( nbQubits , qubits[0] , qubits[1] , control , target ,
             this->controlState() , f ) ;
  }

#ifdef QCLAB_OMP_OFFLOADING
  // apply_device
  template <typename T>
  void CZ< T >::apply_device( Op op , const int nbQubits , T* vector ,
                              const int offset ) const {
    auto qubits = this->qubits() ;
    qubits[0] += offset ;
    qubits[1] += offset ;
    const int control = this->control() + offset ;
    const int target  = this->target()  + offset ;
    auto f = lambda_PauliZ( op , vector ) ;
    apply_device4b( nbQubits , qubits[0] , qubits[1] , control , target ,
                    this->controlState() , f ) ;
  }
#endif

  // apply
  template <typename T>
  void CZ< T >::apply( Side side , Op op , const int nbQubits ,
                       qclab::dense::SquareMatrix< T >& matrix ,
                       const int offset ) const {
    QControlledGate2< T >::apply( side , op , nbQubits , matrix , offset ) ;
  }

  template class CZ< float > ;
  template class CZ< double > ;
  template class CZ< std::complex< float > > ;
  template class CZ< std::complex< double > > ;

} // namespace qclab::qgates

