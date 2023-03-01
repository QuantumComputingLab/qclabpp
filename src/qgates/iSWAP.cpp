#include "qclab/qgates/iSWAP.hpp"
#include "qclab/qgates/Hadamard.hpp"
#include "qclab/qgates/Phase90.hpp"
#include "qclab/qgates/CNOT.hpp"
#include "apply.hpp"

namespace qclab::qgates {

  // apply
  template <typename T>
  void iSWAP< T >::apply( Op op , const int nbQubits ,
                          std::vector< T >& vector , const int offset ) const {
    auto qubits = this->qubits() ;
    qubits[0] += offset ;
    qubits[1] += offset ;
    auto f = lambda_iSWAP( op , vector.data() ) ;
    apply4bc( nbQubits , qubits[0] , qubits[1] , f ) ;
  }

#ifdef QCLAB_OMP_OFFLOADING
  // apply_device
  template <typename T>
  void iSWAP< T >::apply_device( Op op , const int nbQubits , T* vector ,
                                 const int offset ) const {
    auto qubits = this->qubits() ;
    qubits[0] += offset ;
    qubits[1] += offset ;
    auto f = lambda_iSWAP( op , vector ) ;
    apply_device4bc( nbQubits , qubits[0] , qubits[1] , f ) ;
  }
#endif

  template <typename T>
  void iSWAP< T >::apply( Side side , Op op , const int nbQubits ,
                          qclab::dense::SquareMatrix< T >& matrix ,
                          const int offset ) const {
    assert( nbQubits >= 2 ) ;
    assert( matrix.size() == 1 << nbQubits ) ;
    const int qubit0 = qubits_[0] + offset ;
    const int qubit1 = qubits_[1] + offset ;
    assert( qubit0 < nbQubits ) ; assert( qubit1 < nbQubits ) ;
    //
    // -/-----\-   -[S]-[H]--*--[X]-----   -----[X]--*--[H]-[S]-
    //  |iSWAP|  =           |   |       =       |   |
    // -\-----/-   -[S]-----[X]--*--[H]-   -[H]--*--[X]-----[S]-
    //
    qclab::qgates::Phase90< T >  S( qubit0 ) ;
    qclab::qgates::Hadamard< T > H( qubit0 ) ;
    qclab::qgates::CNOT< T >     CNOT( qubit0 , qubit1 ) ;
    // layer 1
    S.apply( side , op , nbQubits , matrix ) ;
    S.setQubit( qubit1 ) ;
    S.apply( side , op , nbQubits , matrix ) ;
    // layer 2
    H.apply( side , op , nbQubits , matrix ) ;
    // layer 3
    CNOT.apply( side , op , nbQubits , matrix ) ;
    // layer 4
    int qnew[] = { qubit1 , qubit0 } ;
    CNOT.setQubits( &qnew[0] ) ;
    CNOT.apply( side , op , nbQubits , matrix ) ;
    // layer 5
    H.setQubit( qubit1 ) ;
    H.apply( side , op , nbQubits , matrix ) ;
  }

  template class iSWAP< std::complex< float > > ;
  template class iSWAP< std::complex< double > > ;

} // namespace qclab::qgates

