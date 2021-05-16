#include "QCircuit.hpp"
#include "qgates/Hadamard.hpp"
#include "qgates/CPhase.hpp"
#include "qgates/SWAP.hpp"
#include <iostream>
#include <iomanip>

template <typename T>
void qft( qclab::QCircuit< T >& circuit ) {

  using R  = qclab::real_t< T > ;
  using H  = qclab::qgates::Hadamard< T > ;
  using CP = qclab::qgates::CPhase< T > ;
  using SW = qclab::qgates::SWAP< T > ;

  // constants
  const R pi = 4 * std::atan(1) ;
  const int n = circuit.nbQubits() ;

  // B blocks
  for ( int i = 0; i < n; i++ ) {
    // Hadamard
    circuit.push_back( std::make_unique< H >( i ) ) ;
    // diagonal blocks
    for ( int j = 2; j <= n-i; j++ ) {
      const int control = j + i - 1 ;
      const R theta = -2*pi / ( 1 << j ) ;
      circuit.push_back( std::make_unique< CP >( control , i , theta ) ) ;
    }
  }

  // swaps
  for ( int i = 0; i < n/2; i++ ) {
    circuit.push_back( std::make_unique< SW >( i , n - i - 1 ) ) ;
  }

}


int main( int argc , char *argv[] ) {

  using T = std::complex< double > ;

  // defaults
  int nbQubits = 2 ;
  int maxPrint = 5 ;

  // arguments
  if ( argc > 1 ) nbQubits = std::stoi( argv[1] ) ;
  if ( argc > 2 ) maxPrint = std::stoi( argv[2] ) ;
  std::cout << "nb qubits = " << nbQubits << std::endl ;

  // quantum circuit
  qclab::QCircuit< T > circuit( nbQubits ) ;

  // qft
  qft( circuit ) ;

  // print matrix
  if ( nbQubits <= maxPrint ) {
    std::cout << "\nmatrix =" << std::endl ;
    qclab::printMatrix( circuit.matrix() ) ;
  }

  // print qasm
  std::stringstream qasm ;
  circuit.toQASM( qasm ) ;
  std::cout << "\nqasm:\n" << qasm.str() ;

  // successful
  return 0 ;

}

