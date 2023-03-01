#include "run.hpp"
#include "qclab/qgates/RotationX.hpp"
#include "qclab/qgates/RotationZ.hpp"
#include "qclab/qgates/CNOT.hpp"

template <typename T>
void trotter( qclab::QCircuit< T >& circuit ) {

  using R  = qclab::real_t< T > ;
  using RX = qclab::qgates::RotationX< T > ;
  using RZ = qclab::qgates::RotationZ< T > ;
  using CX = qclab::qgates::CNOT< T > ;

  // constants
  const R pi = 4 * std::atan(1) ;
  const int n = circuit.nbQubits() ;

  // layers
  for ( int lyr = 0; lyr < n; lyr++ ) {
    int os = ( lyr % 2 ) ? 1 : 0 ;        // offset at bottom
    int os_c = ( n % 2 ) ? 1 - os : os ;  // offset at top
    for ( int i = os; i < n - os_c; i++ ) {
      circuit.push_back( std::make_unique< RZ >( i , pi ) ) ;
    }
    for ( int i = os; i < n - os_c - 1; i += 2 ) {
      circuit.push_back( std::make_unique< RX >( i     , pi/2 ) ) ;
      circuit.push_back( std::make_unique< RX >( i + 1 , pi/2 ) ) ;
      circuit.push_back( std::make_unique< CX >( i , i + 1 ) ) ;
      circuit.push_back( std::make_unique< RX >( i     , pi ) ) ;
      circuit.push_back( std::make_unique< RZ >( i + 1 , pi ) ) ;
      circuit.push_back( std::make_unique< CX >( i , i + 1 ) ) ;
      circuit.push_back( std::make_unique< RX >( i     , -pi/2 ) ) ;
      circuit.push_back( std::make_unique< RX >( i + 1 , -pi/2 ) ) ;
    }
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
    auto f = [&] ( qclab::QCircuit< T >& circuit ) { trotter( circuit ) ; } ;
    r = timings< T >( qmin , qmax , qstp , test , f ) ;
  } else if ( type == 'd' ) {
    // double
    using T = std::complex< double > ;
    auto f = [&] ( qclab::QCircuit< T >& circuit ) { trotter( circuit ) ; } ;
    std::cout << ", T = std::complex<double>" << std::endl ;
    r = timings< T >( qmin , qmax , qstp , test , f ) ;
  } else {
    r = -100 ;
  }
  return r ;

}

