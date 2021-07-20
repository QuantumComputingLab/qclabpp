#include <gtest/gtest.h>
#include "qclab/qgates/Phase45.hpp"

template <typename T>
void test_qclab_qgates_Phase45() {

  using R = qclab::real_t< T > ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  qclab::qgates::Phase45< T >  Tgate ;

  EXPECT_EQ( Tgate.nbQubits() , 1 ) ;   // nbQubits
  EXPECT_TRUE( Tgate.fixed() ) ;        // fixed
  EXPECT_FALSE( Tgate.controlled() ) ;  // controlled

  // qubit
  EXPECT_EQ( Tgate.qubit() , 0 ) ;
  Tgate.setQubit( 2 ) ;
  EXPECT_EQ( Tgate.qubit() , 2 ) ;

  // qubits
  auto qubits = Tgate.qubits() ;
  EXPECT_EQ( qubits.size() , 1 ) ;
  EXPECT_EQ( qubits[0] , 2 ) ;
  int qnew = 3 ;
  Tgate.setQubits( &qnew ) ;
  EXPECT_EQ( Tgate.qubit() , 3 ) ;

  // matrix
  EXPECT_EQ(   std::real( Tgate.matrix()(0,0) ) , 1 ) ;
  EXPECT_EQ(   std::real( Tgate.matrix()(1,0) ) , 0 ) ;
  EXPECT_EQ(   std::real( Tgate.matrix()(0,1) ) , 0 ) ;
  EXPECT_NEAR( std::real( Tgate.matrix()(1,1) ) , 1 / std::sqrt( 2. ) , eps ) ;
  EXPECT_EQ(   std::imag( Tgate.matrix()(0,0) ) , 0 ) ;
  EXPECT_EQ(   std::imag( Tgate.matrix()(1,0) ) , 0 ) ;
  EXPECT_EQ(   std::imag( Tgate.matrix()(0,1) ) , 0 ) ;
  EXPECT_NEAR( std::imag( Tgate.matrix()(1,1) ) , 1 / std::sqrt( 2. ) , eps ) ;

  // print
  Tgate.print() ;

  // toQASM
  std::stringstream qasm ;
  EXPECT_EQ( Tgate.toQASM( qasm ) , 0 ) ;
  EXPECT_EQ( qasm.str() , "t q[3];\n" ) ;
  std::cout << qasm.str() ;

  // operators == and !=
  qclab::qgates::Phase45< T > Tgate2 ;
  EXPECT_TRUE( Tgate == Tgate2 ) ;
  EXPECT_FALSE( Tgate != Tgate2 ) ;

}


/*
 * complex float
 */
TEST( qclab_qgates_Phase45 , complex_float ) {
  test_qclab_qgates_Phase45< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_Phase45 , complex_double ) {
  test_qclab_qgates_Phase45< std::complex< double > >() ;
}

