#include "run.hpp"
#include "qclab/qgates/Hadamard.hpp"
#include "qclab/qgates/CPhase.hpp"
#include "qclab/qgates/SWAP.hpp"

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

  // defaults
  char type = 'd' ;
  int  qmin = 10 ;
  int  qmax = 20 ;
  int  qstp = 2 ;
  int  test = 3 ;

  // arguments
  if ( argc > 1 ) type = argv[1][0] ;
  if ( argc > 2 ) qmin = std::stoi( argv[2] ) ;
  if ( argc > 3 ) qmax = std::stoi( argv[3] ) ;
  if ( argc > 4 ) qstp = std::stoi( argv[4] ) ;
  if ( argc > 5 ) test = std::stoi( argv[5] ) ;
  std::cout << "nb qubits = " << qmin << ":" << qstp << ":" << qmax ;

  int r = 0 ;
  if ( type == 's' ) {
    // float
    std::cout << ", T = std::complex<float>" << std::endl ;
    using T = std::complex< float > ;
    auto f = [&] ( qclab::QCircuit< T >& circuit ) { qft( circuit ) ; } ;
    r = timings< T >( qmin , qmax , qstp , test , f ) ;
  } else if ( type == 'd' ) {
    // double
    using T = std::complex< double > ;
    auto f = [&] ( qclab::QCircuit< T >& circuit ) { qft( circuit ) ; } ;
    std::cout << ", T = std::complex<double>" << std::endl ;
    r = timings< T >( qmin , qmax , qstp , test , f ) ;
  } else {
    r = -100 ;
  }
  return r ;

}

