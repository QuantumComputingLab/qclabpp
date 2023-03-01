#include "qclab/QCircuit.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif

using TP = std::chrono::time_point< std::chrono::high_resolution_clock > ;

inline void tic( TP& t_bgn ) {
  t_bgn = std::chrono::high_resolution_clock::now() ;
}

double toc( TP& t_bgn ) {
  TP t_end = std::chrono::high_resolution_clock::now() ;
  double time = std::chrono::duration< double >(t_end - t_bgn).count() ;
  return time ;
}


template <typename T, typename F>
int run( const int nbQubits , const int test ,
         double& t_cpu , double& t_gpu , F& lambda ,
         const int IMAX_CPU = 3 , const int IMAX_GPU = 3 ) {

  using R = qclab::real_t< T > ;
  const R tol = 10 * std::numeric_limits< R >::epsilon() ;

  // quantum circuit
  qclab::QCircuit< T > circuit( nbQubits ) ;
  lambda( circuit ) ;

  // timing variable
  std::chrono::time_point< std::chrono::high_resolution_clock > time ;
  t_cpu = 9999 ;
  t_gpu = 9999 ;
  int imax_cpu = IMAX_CPU ;
  int imax_gpu = IMAX_GPU ;

  // vector size
  const size_t N = size_t(1) << nbQubits ;

  // simulate on CPU
  std::vector< T > psi ;
  if ( test == 1 || test == 3 ) {
    std::cout << "CPU(" << nbQubits << ")" << std::endl ;
    int i = 0 ;
    while ( i < imax_cpu ) {
      // all zero state
      psi = std::vector< T >( N , T(0) ) ;
      psi[0] = 1 ;
      // simulate
      tic( time ) ;
      circuit.simulate( psi ) ;
      auto t = toc( time ) ;
      if ( t < t_cpu ) t_cpu = t ;
      // increment
      ++i ;
      // update
      if ( i == IMAX_CPU ) {
        if ( t_cpu < 0.1 ) {
          imax_cpu = 100 ;
        } else if ( t_cpu < 1 ) {
          imax_cpu = 10 ;
        }
      }
      // output
      if ( i <= IMAX_CPU || i == imax_cpu ) {
        std::printf( "  * it%4i:%12.6fs" , i , t ) ;
        if ( i != imax_cpu ) std::cout << std::endl ;
      }
    }
    std::printf( "  -->  min t_cpu = %.6fs\n" , t_cpu ) ;
  }
  if ( t_cpu == 9999 ) { t_cpu = 0 ; }

  // simulate on GPU
  if ( test == 2 || test == 3 ) {
  #ifdef QCLAB_OMP_OFFLOADING
    std::cout << "GPU(" << nbQubits << ")" << std::endl ;
    std::vector< T > phi( N ) ; T* phi_ = phi.data() ;
    #pragma omp target enter data map(alloc:phi_[0:N])
    int i = 0 ;
    while ( i < imax_gpu ) {
      // all zero state
      #pragma omp target
      {
        phi_[0] = 1 ;
        #pragma omp parallel for
        for ( size_t i = 1; i < N; i++ ) {
          phi_[i] = 0 ;
        }
      }
      // simulate
      tic( time ) ;
      circuit.simulate_device( phi_ ) ;
      auto t = toc( time ) ;
      if ( t < t_gpu ) t_gpu = t ;
      // increment
      ++i ;
      // update
      if ( i == IMAX_GPU ) {
        if ( t_gpu < 0.1 ) {
          imax_gpu = 100 ;
        } else if ( t_gpu < 1 ) {
          imax_gpu = 10 ;
        }
      }
      // output
      if ( i <= IMAX_GPU || i == imax_gpu ) {
        std::printf( "  * it%4i:%12.6fs" , i , t ) ;
        if ( i != imax_gpu ) std::cout << std::endl ;
      }
    }
    std::printf( "  -->  min t_gpu = %.6fs" , t_gpu ) ;
    // check
    if ( test == 3 ) {
      std::printf( "  -->  speedup = %.2fx\n" , t_cpu / t_gpu ) ;
      #pragma omp target update from(phi_[0:N])
      for ( size_t i = 0; i < N; i++ ) {
        if ( std::abs( psi[i] - phi[i] ) > tol ) {
          std::cout << " *** CHECK FAILED! ***" << std::endl ;
          return -1 ;
        }
      }
      std::cout << "Check passed." << std::endl ;
    } else {
      std::cout << std::endl ;
    }
    #pragma omp target exit data map(delete:phi_[0:N])
  #endif
  }
  if ( t_gpu == 9999 ) { t_gpu = 0 ; }

  // successful
  return 0 ;

}


template <typename T, typename F>
int timings( const int qmin , const int qmax , const int qstep , const int test,
             F& lambda , const int IMAX_CPU = 3 , const int IMAX_GPU = 3 ) {

  // omp
#ifdef _OPENMP
  std::cout << "  --> omp_max_threads() = " << omp_get_max_threads() << "\n" ;
#endif

  // variables
  double t_cpu ;
  double t_gpu ;
  std::vector< int > qubits ;
  std::vector< double > T_cpu ;
  std::vector< double > T_gpu ;

  // run
  for ( int q = qmin; q <= qmax; q += qstep ) {
    int r = run< T , F >( q , test , t_cpu , t_gpu ,
                          lambda , IMAX_CPU , IMAX_GPU ) ;
    qubits.push_back( q ) ;
    T_cpu.push_back( t_cpu ) ;
    T_gpu.push_back( t_gpu ) ;
    if ( r != 0 ) return r ;
  }

  // output
  std::cout << std::endl ;
  for ( int i = 0; i < qubits.size(); ++i ) {
    std::printf( "  %2i | %12.6fs | %12.6fs | %6.2fx\n" ,
                 qubits[i] , T_cpu[i] , T_gpu[i] , T_cpu[i] / T_gpu[i] ) ;
  }
  std::cout << std::endl ;

  // successful
  return 0 ;

}

