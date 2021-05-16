#include <gtest/gtest.h>
#include "qgates/PauliY.hpp"

template <typename T>
void test_qclab_qgates_PauliY() {

  qclab::qgates::PauliY< T >  Y ;

  EXPECT_EQ( Y.nbQubits() , 1 ) ;   // nbQubits
  EXPECT_TRUE( Y.fixed() ) ;        // fixed
  EXPECT_FALSE( Y.controlled() ) ;  // controlled

  // qubit
  EXPECT_EQ( Y.qubit() , 0 ) ;
  Y.setQubit( 2 ) ;
  EXPECT_EQ( Y.qubit() , 2 ) ;

  // qubits
  auto qubits = Y.qubits() ;
  EXPECT_EQ( qubits.size() , 1 ) ;
  EXPECT_EQ( qubits[0] , 2 ) ;
  int qnew = 3 ;
  Y.setQubits( &qnew ) ;
  EXPECT_EQ( Y.qubit() , 3 ) ;

  // matrix
  EXPECT_EQ( std::real( Y.matrix()(0,0) ) ,  0.0 ) ;
  EXPECT_EQ( std::real( Y.matrix()(1,0) ) ,  0.0 ) ;
  EXPECT_EQ( std::real( Y.matrix()(0,1) ) ,  0.0 ) ;
  EXPECT_EQ( std::real( Y.matrix()(1,1) ) ,  0.0 ) ;
  EXPECT_EQ( std::imag( Y.matrix()(0,0) ) ,  0.0 ) ;
  EXPECT_EQ( std::imag( Y.matrix()(1,0) ) ,  1.0 ) ;
  EXPECT_EQ( std::imag( Y.matrix()(0,1) ) , -1.0 ) ;
  EXPECT_EQ( std::imag( Y.matrix()(1,1) ) ,  0.0 ) ;

  // print
  Y.print() ;

  // toQASM
  std::stringstream qasm ;
  EXPECT_EQ( Y.toQASM( qasm ) , 0 ) ;
  EXPECT_EQ( qasm.str() , "y q[3];\n" ) ;
  std::cout << qasm.str() ;

  // operators == and !=
  qclab::qgates::PauliY< T > Y2 ;
  EXPECT_TRUE( Y == Y2 ) ;
  EXPECT_FALSE( Y != Y2 ) ;

}


/*
 * complex float
 */
TEST( qclab_qgates_PauliY , complex_float ) {
  test_qclab_qgates_PauliY< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_PauliY , complex_double ) {
  test_qclab_qgates_PauliY< std::complex< double > >() ;
}

