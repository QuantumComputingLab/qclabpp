#include "qclab/QCircuit.hpp"

int main( int argc , char *argv[] ) {

  using T = std::complex< double > ;

  // qasm file
  if ( argc < 2 ) {
    std::cout << "ERROR: filename is required!" << std::endl ;
    return -1 ;
  }

  // quantum circuit
  qclab::QCircuit< T > circuit( argv[1] ) ;

  // properties
  std::cout << "nbQubits = " << circuit.nbQubits() << std::endl ;
  std::cout << "nbGates  = " << circuit.nbGates() << std::endl ;

  // print qasm
  std::stringstream qasm ;
  circuit.toQASM( qasm ) ;
  std::cout << "\nqasm:\n" << qasm.str() ;

  // successful
  return 0 ;

}

