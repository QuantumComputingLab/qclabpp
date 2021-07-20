#include <gtest/gtest.h>
#include "qclab/qgates/PauliZ.hpp"

template <typename T>
void test_qclab_qgates_PauliZ() {

  qclab::qgates::PauliZ< T >  Z ;

  EXPECT_EQ( Z.nbQubits() , 1 ) ;   // nbQubits
  EXPECT_TRUE( Z.fixed() ) ;        // fixed
  EXPECT_FALSE( Z.controlled() ) ;  // controlled

  // qubit
  EXPECT_EQ( Z.qubit() , 0 ) ;
  Z.setQubit( 2 ) ;
  EXPECT_EQ( Z.qubit() , 2 ) ;

  // qubits
  auto qubits = Z.qubits() ;
  EXPECT_EQ( qubits.size() , 1 ) ;
  EXPECT_EQ( qubits[0] , 2 ) ;
  int qnew = 3 ;
  Z.setQubits( &qnew ) ;
  EXPECT_EQ( Z.qubit() , 3 ) ;

  // matrix
  EXPECT_EQ( Z.matrix()(0,0) , T( 1.0) ) ;
  EXPECT_EQ( Z.matrix()(1,0) , T( 0.0) ) ;
  EXPECT_EQ( Z.matrix()(0,1) , T( 0.0) ) ;
  EXPECT_EQ( Z.matrix()(1,1) , T(-1.0) ) ;

  // print
  Z.print() ;

  // toQASM
  std::stringstream qasm ;
  EXPECT_EQ( Z.toQASM( qasm ) , 0 ) ;
  EXPECT_EQ( qasm.str() , "z q[3];\n" ) ;
  std::cout << qasm.str() ;

  // operators == and !=
  qclab::qgates::PauliZ< T > Z2 ;
  EXPECT_TRUE( Z == Z2 ) ;
  EXPECT_FALSE( Z != Z2 ) ;

}


/*
 * float
 */
TEST( qclab_qgates_PauliZ , float ) {
  test_qclab_qgates_PauliZ< float >() ;
}

/*
 * double
 */
TEST( qclab_qgates_PauliZ , double ) {
  test_qclab_qgates_PauliZ< double >() ;
}


/*
 * complex float
 */
TEST( qclab_qgates_PauliZ , complex_float ) {
  test_qclab_qgates_PauliZ< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_PauliZ , complex_double ) {
  test_qclab_qgates_PauliZ< std::complex< double > >() ;
}

