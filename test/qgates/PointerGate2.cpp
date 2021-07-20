#include <gtest/gtest.h>
#include "qclab/qgates/PointerGate2.hpp"
#include "qclab/qgates/CNOT.hpp"
#include "qclab/qgates/SWAP.hpp"

template <typename T>
void test_qclab_qgates_PointerGate2() {

  {
    qclab::qgates::SWAP< T >  swap ;
    qclab::qgates::PointerGate2< T >  P( &swap )  ;

    EXPECT_EQ( P.nbQubits() , 2 ) ;   // nbQubits
    EXPECT_TRUE( P.fixed() ) ;        // fixed
    EXPECT_FALSE( P.controlled() ) ;  // controlled
    EXPECT_EQ( P.qubit() , 0 ) ;      // qubit
    EXPECT_EQ( P.offset() , 0 ) ;     // offset
    EXPECT_EQ( P.ptr() , &swap ) ;    // ptr

    // qubits
    auto qubits = P.qubits() ;
    EXPECT_EQ( qubits.size() , 2 ) ;
    EXPECT_EQ( qubits[0] , 0 ) ;
    EXPECT_EQ( qubits[1] , 1 ) ;

    // matrix
    EXPECT_TRUE( P.matrix() == swap.matrix() ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( P.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "swap q[0], q[1];\n" ) ;

    // setOffset
    P.setOffset( 3 ) ;
    EXPECT_EQ( P.offset() , 3 ) ;
    EXPECT_EQ( P.qubits()[0] , 3 ) ;
    EXPECT_EQ( P.qubits()[1] , 4 ) ;
    std::stringstream qasm2 ;
    EXPECT_EQ( P.toQASM( qasm2 ) , 0 ) ;
    EXPECT_EQ( qasm2.str() , "swap q[3], q[4];\n" ) ;

    // operators == and !=
    EXPECT_TRUE(  P == swap ) ;
    EXPECT_FALSE( P != swap ) ;
    qclab::qgates::CNOT< T >  cnot ;
    EXPECT_TRUE(  P != cnot ) ;
    EXPECT_FALSE( P == cnot ) ;
    qclab::qgates::SWAP< T >  swap2 ;
    qclab::qgates::PointerGate2< T > Pswap( &swap2 ) ;
    qclab::qgates::PointerGate2< T > Pcnot( &cnot ) ;
    EXPECT_TRUE(  P == Pswap ) ;
    EXPECT_FALSE( P != Pswap ) ;
    EXPECT_TRUE(  P != Pcnot ) ;
    EXPECT_FALSE( P == Pcnot ) ;
  }

  {
    qclab::qgates::SWAP< T >  swap( 2 , 5 ) ;
    qclab::qgates::PointerGate2< T >  P( &swap )  ;

    EXPECT_EQ( P.nbQubits() , 2 ) ;   // nbQubits
    EXPECT_TRUE( P.fixed() ) ;        // fixed
    EXPECT_FALSE( P.controlled() ) ;  // controlled
    EXPECT_EQ( P.qubits()[0] , 2 ) ;  // qubit0
    EXPECT_EQ( P.qubits()[1] , 5 ) ;  // qubit1
    EXPECT_EQ( P.offset() , 0 ) ;     // offset
    EXPECT_EQ( P.ptr() , &swap ) ;    // ptr

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( P.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "swap q[2], q[5];\n" ) ;
  }

  {
    qclab::qgates::SWAP< T >  swap ;
    qclab::qgates::PointerGate2< T >  P( &swap , 3 )  ;

    EXPECT_EQ( P.nbQubits() , 2 ) ;   // nbQubits
    EXPECT_TRUE( P.fixed() ) ;        // fixed
    EXPECT_FALSE( P.controlled() ) ;  // controlled
    EXPECT_EQ( P.qubits()[0] , 3 ) ;  // qubit0
    EXPECT_EQ( P.qubits()[1] , 4 ) ;  // qubit1
    EXPECT_EQ( P.offset() , 3 ) ;     // offset
    EXPECT_EQ( P.ptr() , &swap ) ;    // ptr

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( P.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "swap q[3], q[4];\n" ) ;
  }

  {
    qclab::qgates::SWAP< T >  swap( 2 , 5 ) ;
    qclab::qgates::PointerGate2< T >  P( &swap , 2 )  ;

    EXPECT_EQ( P.nbQubits() , 2 ) ;   // nbQubits
    EXPECT_TRUE( P.fixed() ) ;        // fixed
    EXPECT_FALSE( P.controlled() ) ;  // controlled
    EXPECT_EQ( P.qubits()[0] , 4 ) ;  // qubit0
    EXPECT_EQ( P.qubits()[1] , 7 ) ;  // qubit1
    EXPECT_EQ( P.offset() , 2 ) ;     // offset
    EXPECT_EQ( P.ptr() , &swap ) ;    // ptr

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( P.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "swap q[4], q[7];\n" ) ;
  }

}


/*
 * float
 */
TEST( qclab_qgates_PointerGate2 , float ) {
  test_qclab_qgates_PointerGate2< float >() ;
}

/*
 * double
 */
TEST( qclab_qgates_PointerGate2 , double ) {
  test_qclab_qgates_PointerGate2< double >() ;
}


/*
 * complex float
 */
TEST( qclab_qgates_PointerGate2 , complex_float ) {
  test_qclab_qgates_PointerGate2< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_PointerGate2 , complex_double ) {
  test_qclab_qgates_PointerGate2< std::complex< double > >() ;
}

