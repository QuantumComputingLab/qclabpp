#include <gtest/gtest.h>
#include "qclab/qgates/PauliX.hpp"

template <typename T>
void test_qclab_qgates_PauliX() {

  qclab::qgates::PauliX< T >  X ;

  EXPECT_EQ( X.nbQubits() , 1 ) ;   // nbQubits
  EXPECT_TRUE( X.fixed() ) ;        // fixed
  EXPECT_FALSE( X.controlled() ) ;  // controlled

  // qubit
  EXPECT_EQ( X.qubit() , 0 ) ;
  X.setQubit( 2 ) ;
  EXPECT_EQ( X.qubit() , 2 ) ;

  // qubits
  auto qubits = X.qubits() ;
  EXPECT_EQ( qubits.size() , 1 ) ;
  EXPECT_EQ( qubits[0] , 2 ) ;
  int qnew = 3 ;
  X.setQubits( &qnew ) ;
  EXPECT_EQ( X.qubit() , 3 ) ;

  // matrix
  EXPECT_EQ( X.matrix()(0,0) , T(0.0) ) ;
  EXPECT_EQ( X.matrix()(1,0) , T(1.0) ) ;
  EXPECT_EQ( X.matrix()(0,1) , T(1.0) ) ;
  EXPECT_EQ( X.matrix()(1,1) , T(0.0) ) ;

  // print
  X.print() ;

  // toQASM
  std::stringstream qasm ;
  EXPECT_EQ( X.toQASM( qasm ) , 0 ) ;
  EXPECT_EQ( qasm.str() , "x q[3];\n" ) ;
  std::cout << qasm.str() ;

  // operators == and !=
  qclab::qgates::PauliX< T > X2 ;
  EXPECT_TRUE( X == X2 ) ;
  EXPECT_FALSE( X != X2 ) ;

}


/*
 * float
 */
TEST( qclab_qgates_PauliX , float ) {
  test_qclab_qgates_PauliX< float >() ;
}

/*
 * double
 */
TEST( qclab_qgates_PauliX , double ) {
  test_qclab_qgates_PauliX< double >() ;
}


/*
 * complex float
 */
TEST( qclab_qgates_PauliX , complex_float ) {
  test_qclab_qgates_PauliX< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_PauliX , complex_double ) {
  test_qclab_qgates_PauliX< std::complex< double > >() ;
}

