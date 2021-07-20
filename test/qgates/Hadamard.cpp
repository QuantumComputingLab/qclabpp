#include <gtest/gtest.h>
#include "qclab/qgates/Hadamard.hpp"
#include <limits>

template <typename T>
void test_qclab_qgates_Hadamard() {

  using R = qclab::real_t< T > ;
  const R eps = std::numeric_limits< R >::epsilon() ;

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
  EXPECT_NEAR( std::real( H.matrix()(0,0) ) ,  1.0 / std::sqrt( 2.0 ) , eps ) ;
  EXPECT_NEAR( std::real( H.matrix()(1,0) ) ,  1.0 / std::sqrt( 2.0 ) , eps ) ;
  EXPECT_NEAR( std::real( H.matrix()(0,1) ) ,  1.0 / std::sqrt( 2.0 ) , eps ) ;
  EXPECT_NEAR( std::real( H.matrix()(1,1) ) , -1.0 / std::sqrt( 2.0 ) , eps ) ;
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

