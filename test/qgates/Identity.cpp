#include <gtest/gtest.h>
#include "qclab/qgates/Identity.hpp"

template <typename T>
void test_qclab_qgates_Identity() {

  qclab::qgates::Identity< T >  I ;

  EXPECT_EQ( I.nbQubits() , 1 ) ;   // nbQubits
  EXPECT_TRUE( I.fixed() ) ;        // fixed
  EXPECT_FALSE( I.controlled() ) ;  // controlled

  // qubit
  EXPECT_EQ( I.qubit() , 0 ) ;
  I.setQubit( 2 ) ;
  EXPECT_EQ( I.qubit() , 2 ) ;

  // qubits
  auto qubits = I.qubits() ;
  EXPECT_EQ( qubits.size() , 1 ) ;
  EXPECT_EQ( qubits[0] , 2 ) ;
  int qnew = 3 ;
  I.setQubits( &qnew ) ;
  EXPECT_EQ( I.qubit() , 3 ) ;

  // matrix
  EXPECT_EQ( I.matrix()(0,0) , T(1.0) ) ;
  EXPECT_EQ( I.matrix()(1,0) , T(0.0) ) ;
  EXPECT_EQ( I.matrix()(0,1) , T(0.0) ) ;
  EXPECT_EQ( I.matrix()(1,1) , T(1.0) ) ;

  // print
  I.print() ;

  // toQASM
  std::stringstream qasm ;
  EXPECT_EQ( I.toQASM( qasm ) , 0 ) ;
  EXPECT_EQ( qasm.str() , "" ) ;

  // operators == and !=
  qclab::qgates::Identity< T > I2 ;
  EXPECT_TRUE( I == I2 ) ;
  EXPECT_FALSE( I != I2 ) ;

}


/*
 * float
 */
TEST( qclab_qgates_Identity , float ) {
  test_qclab_qgates_Identity< float >() ;
}

/*
 * double
 */
TEST( qclab_qgates_Identity , double ) {
  test_qclab_qgates_Identity< double >() ;
}


/*
 * complex float
 */
TEST( qclab_qgates_Identity , complex_float ) {
  test_qclab_qgates_Identity< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_Identity , complex_double ) {
  test_qclab_qgates_Identity< std::complex< double > >() ;
}

