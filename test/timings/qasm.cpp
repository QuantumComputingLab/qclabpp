#include "qclab/QCircuit.hpp"
#include "qclab/qgates/Hadamard.hpp"
#include "qclab/qgates/PauliX.hpp"
#include "qclab/qgates/PauliY.hpp"
#include "qclab/qgates/PauliZ.hpp"
#include "qclab/qgates/RotationX.hpp"
#include "qclab/qgates/RotationY.hpp"
#include "qclab/qgates/RotationZ.hpp"
#include "qclab/qgates/Phase.hpp"
#include "qclab/qgates/Phase45.hpp"
#include "qclab/qgates/Phase90.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>

using TP = std::chrono::time_point< std::chrono::high_resolution_clock > ;

void tic( const std::string name , TP& t_bgn ) {
  std::cout << name << std::endl ;
  t_bgn = std::chrono::high_resolution_clock::now() ;
}

void toc( const TP& t_bgn ) {
  TP t_end = std::chrono::high_resolution_clock::now() ;
  double time = std::chrono::duration< double >(t_end - t_bgn).count() ;
  std::cout << "                          " << time << "s" << std::endl ;
}


template <typename T>
int timings( const std::string& filename , const int test ) {

  using R = qclab::real_t< T > ;
  const R tol = 10 * std::numeric_limits< R >::epsilon() ;

  // timing variable
  std::chrono::time_point< std::chrono::high_resolution_clock > time ;

  // quantum circuit
  qclab::QCircuit< T > circuit( filename ) ;
  std::stringstream qasm ;
  circuit.toQASM( qasm ) ;
  std::cout << "\nqasm:\n" << qasm.str() << std::endl ;

  // vector size
  const size_t N = size_t(1) << circuit.nbQubits() ;

  // simulate on CPU
  std::vector< T > psi ;
  if ( test == 1 || test == 3 ) {
    std::cout << "CPU" << std::endl ;
    // all zero state
    tic( "  * initialize ket0:" , time ) ;
    psi = std::vector< T >( N , T(0) ) ;
    psi[0] = 1 ;
    toc( time ) ;
    // simulate
    tic( "  * simulate qft circuit:" , time ) ;
    circuit.simulate( psi ) ;
    toc( time ) ;
  }

  // simulate on GPU
  if ( test == 2 || test == 3 ) {
  #ifdef QCLAB_OMP_OFFLOADING
    std::cout << "GPU" << std::endl ;
    std::vector< T > phi( N ) ; T* phi_ = phi.data() ;
    #pragma omp target enter data map(alloc:phi_[0:N])
    // all zero state
    tic( "  * initialize ket0:" , time ) ;
    #pragma omp target
    {
      phi_[0] = 1 ;
      #pragma omp parallel for
      for ( size_t i = 1; i < N; i++ ) {
        phi_[i] = 0 ;
      }
    }
    toc( time ) ;
    // simulate
    tic( "  * simulate qft circuit:" , time ) ;
    circuit.simulate_device( phi_ ) ;
    toc( time ) ;
    // check
    if ( test == 3 ) {
      #pragma omp target update from(phi_[0:N])
      #pragma omp target exit data map(delete:phi_[0:N])
      for ( size_t i = 0; i < N; i++ ) {
        if ( std::abs( psi[i] - phi[i] ) > tol ) {
          std::cout << " *** CHECK FAILED! ***" << std::endl ;
          return -1 ;
        }
      }
      std::cout << "Check passed." << std::endl ;
    }
  #endif
  }

  // successful
  return 0 ;

}


int main( int argc , char *argv[] ) {

  // qasm file
  if ( argc < 2 ) {
    std::cout << "ERROR: filename is required!" << std::endl ;
    return -1 ;
  }
  std::string filename( argv[1] ) ;

  // defaults
  char type = 'd' ;
  int  test = 3 ;

  // arguments
  if ( argc > 2 ) type = argv[2][0] ;
  if ( argc > 3 ) test = std::stoi( argv[3] ) ;

  int r = 0 ;
  if ( type == 's' ) {
    // float
    std::cout << "T = std::complex<float>" << std::endl ;
    r = timings< std::complex< float > >( filename , test ) ;
  } else if ( type == 'd' ) {
    // double
    std::cout << "T = std::complex<double>" << std::endl ;
    r = timings< std::complex< double > >( filename , test ) ;
  } else {
    r = -100 ;
  }
  return r ;

}

