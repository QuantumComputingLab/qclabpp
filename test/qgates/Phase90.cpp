#include <gtest/gtest.h>
#include "qclab/qgates/Phase90.hpp"

template <typename T>
void test_qclab_qgates_Phase90() {

  qclab::qgates::Phase90< T >  S ;

  EXPECT_EQ( S.nbQubits() , 1 ) ;   // nbQubits
  EXPECT_TRUE( S.fixed() ) ;        // fixed
  EXPECT_FALSE( S.controlled() ) ;  // controlled

  // qubit
  EXPECT_EQ( S.qubit() , 0 ) ;
  S.setQubit( 2 ) ;
  EXPECT_EQ( S.qubit() , 2 ) ;

  // qubits
  auto qubits = S.qubits() ;
  EXPECT_EQ( qubits.size() , 1 ) ;
  EXPECT_EQ( qubits[0] , 2 ) ;
  int qnew = 3 ;
  S.setQubits( &qnew ) ;
  EXPECT_EQ( S.qubit() , 3 ) ;

  // matrix
  EXPECT_EQ( std::real( S.matrix()(0,0) ) , 1.0 ) ;
  EXPECT_EQ( std::real( S.matrix()(1,0) ) , 0.0 ) ;
  EXPECT_EQ( std::real( S.matrix()(0,1) ) , 0.0 ) ;
  EXPECT_EQ( std::real( S.matrix()(1,1) ) , 0.0 ) ;
  EXPECT_EQ( std::imag( S.matrix()(0,0) ) , 0.0 ) ;
  EXPECT_EQ( std::imag( S.matrix()(1,0) ) , 0.0 ) ;
  EXPECT_EQ( std::imag( S.matrix()(0,1) ) , 0.0 ) ;
  EXPECT_EQ( std::imag( S.matrix()(1,1) ) , 1.0 ) ;

  // print
  S.print() ;

  // toQASM
  std::stringstream qasm ;
  EXPECT_EQ( S.toQASM( qasm ) , 0 ) ;
  EXPECT_EQ( qasm.str() , "s q[3];\n" ) ;
  std::cout << qasm.str() ;

  // operators == and !=
  qclab::qgates::Phase90< T > S2 ;
  EXPECT_TRUE( S == S2 ) ;
  EXPECT_FALSE( S != S2 ) ;

}


/*
 * complex float
 */
TEST( qclab_qgates_Phase90 , complex_float ) {
  test_qclab_qgates_Phase90< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_Phase90 , complex_double ) {
  test_qclab_qgates_Phase90< std::complex< double > >() ;
}

