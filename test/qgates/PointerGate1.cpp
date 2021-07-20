#include <gtest/gtest.h>
#include "qclab/qgates/PointerGate1.hpp"
#include "qclab/qgates/PauliX.hpp"
#include "qclab/qgates/PauliZ.hpp"

template <typename T>
void test_qclab_qgates_PointerGate1() {

  {
    qclab::qgates::PauliX< T >  X ;
    qclab::qgates::PointerGate1< T >  P( &X )  ;

    EXPECT_EQ( P.nbQubits() , 1 ) ;   // nbQubits
    EXPECT_TRUE( P.fixed() ) ;        // fixed
    EXPECT_FALSE( P.controlled() ) ;  // controlled
    EXPECT_EQ( P.qubit() , 0 ) ;      // qubit
    EXPECT_EQ( P.offset() , 0 ) ;     // offset
    EXPECT_EQ( P.ptr() , &X ) ;       // ptr

    // qubits
    auto qubits = P.qubits() ;
    EXPECT_EQ( qubits.size() , 1 ) ;
    EXPECT_EQ( qubits[0] , 0 ) ;

    // matrix
    EXPECT_TRUE( P.matrix() == X.matrix() ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( P.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[0];\n" ) ;

    // setOffset
    P.setOffset( 3 ) ;
    EXPECT_EQ( P.offset() , 3 ) ;
    EXPECT_EQ( P.qubit() , 3 ) ;
    EXPECT_EQ( P.qubits()[0] , 3 ) ;
    std::stringstream qasm2 ;
    EXPECT_EQ( P.toQASM( qasm2 ) , 0 ) ;
    EXPECT_EQ( qasm2.str() , "x q[3];\n" ) ;

    // operators == and !=
    EXPECT_TRUE(  P == X ) ;
    EXPECT_FALSE( P != X ) ;
    qclab::qgates::PauliZ< T >  Z ;
    EXPECT_TRUE(  P != Z ) ;
    EXPECT_FALSE( P == Z ) ;
    qclab::qgates::PauliX< T >  X2 ;
    qclab::qgates::PointerGate1< T > PX( &X2 ) ;
    qclab::qgates::PointerGate1< T > PZ( &Z ) ;
    EXPECT_TRUE(  P == PX ) ;
    EXPECT_FALSE( P != PX ) ;
    EXPECT_TRUE(  P != PZ ) ;
    EXPECT_FALSE( P == PZ ) ;
  }

  {
    qclab::qgates::PauliX< T >  X( 5 ) ;
    qclab::qgates::PointerGate1< T >  P( &X )  ;

    EXPECT_EQ( P.nbQubits() , 1 ) ;   // nbQubits
    EXPECT_TRUE( P.fixed() ) ;        // fixed
    EXPECT_FALSE( P.controlled() ) ;  // controlled
    EXPECT_EQ( P.qubit() , 5 ) ;      // qubit
    EXPECT_EQ( P.qubits()[0] , 5 ) ;  // qubits
    EXPECT_EQ( P.offset() , 0 ) ;     // offset
    EXPECT_EQ( P.ptr() , &X ) ;       // ptr

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( P.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[5];\n" ) ;
  }

  {
    qclab::qgates::PauliX< T >  X ;
    qclab::qgates::PointerGate1< T >  P( &X , 3 )  ;

    EXPECT_EQ( P.nbQubits() , 1 ) ;   // nbQubits
    EXPECT_TRUE( P.fixed() ) ;        // fixed
    EXPECT_FALSE( P.controlled() ) ;  // controlled
    EXPECT_EQ( P.qubit() , 3 ) ;      // qubit
    EXPECT_EQ( P.qubits()[0] , 3 ) ;  // qubits
    EXPECT_EQ( P.offset() , 3 ) ;     // offset
    EXPECT_EQ( P.ptr() , &X ) ;       // ptr

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( P.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[3];\n" ) ;
  }

  {
    qclab::qgates::PauliX< T >  X( 4 ) ;
    qclab::qgates::PointerGate1< T >  P( &X , 3 )  ;

    EXPECT_EQ( P.nbQubits() , 1 ) ;   // nbQubits
    EXPECT_TRUE( P.fixed() ) ;        // fixed
    EXPECT_FALSE( P.controlled() ) ;  // controlled
    EXPECT_EQ( P.qubit() , 7 ) ;      // qubit
    EXPECT_EQ( P.qubits()[0] , 7 ) ;  // qubits
    EXPECT_EQ( P.offset() , 3 ) ;     // offset
    EXPECT_EQ( P.ptr() , &X ) ;       // ptr

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( P.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[7];\n" ) ;
  }

}


/*
 * float
 */
TEST( qclab_qgates_PointerGate1 , float ) {
  test_qclab_qgates_PointerGate1< float >() ;
}

/*
 * double
 */
TEST( qclab_qgates_PointerGate1 , double ) {
  test_qclab_qgates_PointerGate1< double >() ;
}


/*
 * complex float
 */
TEST( qclab_qgates_PointerGate1 , complex_float ) {
  test_qclab_qgates_PointerGate1< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_PointerGate1 , complex_double ) {
  test_qclab_qgates_PointerGate1< std::complex< double > >() ;
}

