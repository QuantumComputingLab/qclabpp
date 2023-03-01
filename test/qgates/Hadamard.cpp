#include <gtest/gtest.h>
#include "qclab/qgates/Hadamard.hpp"
#include <limits>

template <typename T>
void check_Hadamard( const std::vector< T >& v1 , const std::vector< T >& v2 ) {

  using R = qclab::real_t< T > ;
  const R tol = 10 * std::numeric_limits< R >::epsilon() ;

  assert( v1.size() == v2.size() ) ;
  for ( size_t i = 0; i < v1.size(); ++i ) {
    EXPECT_NEAR( std::real( v1[i] ) , std::real( v2[i] ) , tol ) ;
    EXPECT_NEAR( std::imag( v1[i] ) , std::imag( v2[i] ) , tol ) ;
  }

}

template <typename T>
void test_qclab_qgates_Hadamard() {

  using R = qclab::real_t< T > ;
  const R tol = 10 * std::numeric_limits< R >::epsilon() ;

  qclab::qgates::Hadamard< T >  H ;

  EXPECT_EQ( H.nbQubits() , 1 ) ;   // nbQubits
  EXPECT_TRUE( H.fixed() ) ;        // fixed
  EXPECT_FALSE( H.controlled() ) ;  // controlled

  // qubit
  EXPECT_EQ( H.qubit() , 0 ) ;
  H.setQubit( 2 ) ;
  EXPECT_EQ( H.qubit() , 2 ) ;

  // qubits
  auto qubits = H.qubits() ;
  EXPECT_EQ( qubits.size() , 1 ) ;
  EXPECT_EQ( qubits[0] , 2 ) ;
  int qnew = 3 ;
  H.setQubits( &qnew ) ;
  EXPECT_EQ( H.qubit() , 3 ) ;

  // matrix
  EXPECT_NEAR( std::real( H.matrix()(0,0) ) ,  1.0 / std::sqrt( 2.0 ) , tol ) ;
  EXPECT_NEAR( std::real( H.matrix()(1,0) ) ,  1.0 / std::sqrt( 2.0 ) , tol ) ;
  EXPECT_NEAR( std::real( H.matrix()(0,1) ) ,  1.0 / std::sqrt( 2.0 ) , tol ) ;
  EXPECT_NEAR( std::real( H.matrix()(1,1) ) , -1.0 / std::sqrt( 2.0 ) , tol ) ;
  EXPECT_EQ( std::imag( H.matrix()(0,0) ) , 0.0 ) ;
  EXPECT_EQ( std::imag( H.matrix()(1,0) ) , 0.0 ) ;
  EXPECT_EQ( std::imag( H.matrix()(0,1) ) , 0.0 ) ;
  EXPECT_EQ( std::imag( H.matrix()(1,1) ) , 0.0 ) ;

  // print
  H.print() ;

  // toQASM
  std::stringstream qasm ;
  EXPECT_EQ( H.toQASM( qasm ) , 0 ) ;
  EXPECT_EQ( qasm.str() , "h q[3];\n" ) ;
  std::cout << qasm.str() ;

  // operators == and !=
  qclab::qgates::Hadamard< T > H2 ;
  EXPECT_TRUE( H == H2 ) ;
  EXPECT_FALSE( H != H2 ) ;

  // apply
  {
    const R c = R(1) / std::sqrt( R(2) ) ;
    using V = std::vector< T > ;
    const V v1 = { 3 , 5 } ;
    const V v2 = { 3 , 5 , 2 , 7 } ;

    qclab::qgates::Hadamard< T >  H0( 0 ) ;
    qclab::qgates::Hadamard< T >  H1( 1 ) ;

    // nbQubits = 1
    auto vec1 = v1 ;
    H0.apply( qclab::Op::NoTrans , 1 , vec1 ) ;
    V check1 = { c*T(3) + c*T(5) , c*T(3) - c*T(5) } ;
    check_Hadamard( vec1 , check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ; T* vec1_ = vec1.data() ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { H0.apply_device( qclab::Op::NoTrans , 1 , vec1_ ) ; }
    check_Hadamard( vec1 , check1 ) ;
  #endif

    vec1 = v1 ;
    H0.apply( qclab::Op::Trans , 1 , vec1 ) ;
    check_Hadamard( vec1 , check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { H0.apply_device( qclab::Op::Trans , 1 , vec1_ ) ; }
    check_Hadamard( vec1 , check1 ) ;
  #endif

    vec1 = v1 ;
    H0.apply( qclab::Op::ConjTrans , 1 , vec1 ) ;
    check_Hadamard( vec1 , check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { H0.apply_device( qclab::Op::ConjTrans , 1 , vec1_ ) ; }
    check_Hadamard( vec1 , check1 ) ;
  #endif

    // nbQubits = 2
    auto vec2 = v2 ;
    H0.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { c*T(3) + c*T(2) , c*T(5) + c*T(7) ,
                 c*T(3) - c*T(2) , c*T(5) - c*T(7) } ;
    check_Hadamard( vec2 , check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { H0.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    check_Hadamard( vec2 , check2 ) ;
  #endif

    vec2 = v2 ;
    H0.apply( qclab::Op::Trans , 2 , vec2 ) ;
    check_Hadamard( vec2 , check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { H0.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    check_Hadamard( vec2 , check2 ) ;
  #endif

    vec2 = v2 ;
    H0.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check_Hadamard( vec2 , check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { H0.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    check_Hadamard( vec2 , check2 ) ;
  #endif

    vec2 = v2 ;
    H1.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    check2 = { c*T(3) + c*T(5) , c*T(3) - c*T(5) ,
               c*T(2) + c*T(7) , c*T(2) - c*T(7) } ;
    check_Hadamard( vec2 , check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { H1.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    check_Hadamard( vec2 , check2 ) ;
  #endif

    vec2 = v2 ;
    H1.apply( qclab::Op::Trans , 2 , vec2 ) ;
    check_Hadamard( vec2 , check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { H1.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    check_Hadamard( vec2 , check2 ) ;
  #endif

    vec2 = v2 ;
    H1.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check_Hadamard( vec2 , check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { H1.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    check_Hadamard( vec2 , check2 ) ;
  #endif
  }

}


/*
 * float
 */
TEST( qclab_qgates_Hadamard , float ) {
  test_qclab_qgates_Hadamard< float >() ;
}

/*
 * double
 */
TEST( qclab_qgates_Hadamard , double ) {
  test_qclab_qgates_Hadamard< double >() ;
}


/*
 * complex float
 */
TEST( qclab_qgates_Hadamard , complex_float ) {
  test_qclab_qgates_Hadamard< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_Hadamard , complex_double ) {
  test_qclab_qgates_Hadamard< std::complex< double > >() ;
}

