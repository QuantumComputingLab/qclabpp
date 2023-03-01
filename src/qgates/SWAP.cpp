#include "qclab/qgates/SWAP.hpp"
#include "qclab/qgates/CNOT.hpp"
#include "apply.hpp"

namespace qclab::qgates {

  // apply
  template <typename T>
  void SWAP< T >::apply( Op op , const int nbQubits , std::vector< T >& vector ,
                         const int offset ) const {
    auto qubits = this->qubits() ;
    qubits[0] += offset ;
    qubits[1] += offset ;
    auto f = lambda_SWAP( op , vector.data() ) ;
    apply4bc( nbQubits , qubits[0] , qubits[1] , f ) ;
  }

#ifdef QCLAB_OMP_OFFLOADING
  // apply_device
  template <typename T>
  void SWAP< T >::apply_device( Op op , const int nbQubits , T* vector ,
                                const int offset ) const {
    auto qubits = this->qubits() ;
    qubits[0] += offset ;
    qubits[1] += offset ;
    auto f = lambda_SWAP( op , vector ) ;
    apply_device4bc( nbQubits , qubits[0] , qubits[1] , f ) ;
  }
#endif

  // apply
  template <typename T>
  void SWAP< T >::apply( Side side , Op op , const int nbQubits ,
                         qclab::dense::SquareMatrix< T >& matrix ,
                         const int offset ) const {
    assert( nbQubits >= 2 ) ;
    assert( matrix.size() == 1 << nbQubits ) ;
    const int qubit0 = qubits_[0] + offset ;
    const int qubit1 = qubits_[1] + offset ;
    assert( qubit0 < nbQubits ) ; assert( qubit1 < nbQubits ) ;
    // 3x CNOT
    qclab::qgates::CNOT< T > cnot01( qubit0 , qubit1 ) ;
    qclab::qgates::CNOT< T > cnot10( qubit1 , qubit0 ) ;
    cnot01.apply( side , op , nbQubits , matrix ) ;
    cnot10.apply( side , op , nbQubits , matrix ) ;
    cnot01.apply( side , op , nbQubits , matrix ) ;
  }

  template class SWAP< float > ;
  template class SWAP< double > ;
  template class SWAP< std::complex< float > > ;
  template class SWAP< std::complex< double > > ;

} // namespace qclab::qgates

