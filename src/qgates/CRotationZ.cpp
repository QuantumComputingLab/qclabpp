#include "qclab/qgates/CRotationZ.hpp"
#include "apply.hpp"

namespace qclab::qgates {

  // apply
  template <typename T>
  void CRotationZ< T >::apply( Op op , const int nbQubits ,
                               std::vector< T >& vector ,
                               const int offset ) const {
    auto qubits = this->qubits() ;
    qubits[0] += offset ;
    qubits[1] += offset ;
    const int control = this->control() + offset ;
    const int target  = this->target()  + offset ;
    auto f = lambda_RotationZ( op , this->cos() , this->sin() , vector.data() );
    apply4( nbQubits , qubits[0] , qubits[1] , control , target ,
            this->controlState() , f ) ;
  }

#ifdef QCLAB_OMP_OFFLOADING
  // apply_device
  template <typename T>
  void CRotationZ< T >::apply_device( Op op , const int nbQubits , T* vector ,
                                      const int offset ) const {
    auto qubits = this->qubits() ;
    qubits[0] += offset ;
    qubits[1] += offset ;
    const int control = this->control() + offset ;
    const int target  = this->target()  + offset ;
    auto f = lambda_RotationZ( op , this->cos() , this->sin() , vector ) ;
    apply_device4( nbQubits , qubits[0] , qubits[1] , control , target ,
                   this->controlState() , f ) ;
  }
#endif

  // apply
  template <typename T>
  void CRotationZ< T >::apply( Side side , Op op , const int nbQubits ,
                               qclab::dense::SquareMatrix< T >& matrix ,
                               const int offset ) const {
    QControlledGate2< T >::apply( side , op , nbQubits , matrix , offset ) ;
  }

  template class CRotationZ< std::complex< float > > ;
  template class CRotationZ< std::complex< double > > ;

} // namespace qclab::qgates

