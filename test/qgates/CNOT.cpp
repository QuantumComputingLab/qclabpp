#include <gtest/gtest.h>
#include "qclab/qgates/CNOT.hpp"

template <typename T>
void test_qclab_qgates_CNOT() {

  //
  // CNOTs |1> controlled
  //
  {
    qclab::qgates::CNOT< T >  cnot ;

    EXPECT_EQ( cnot.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cnot.fixed() ) ;           // fixed
    EXPECT_TRUE( cnot.controlled() ) ;      // controlled
    EXPECT_EQ( cnot.control() , 0 ) ;       // control
    EXPECT_EQ( cnot.target() , 1 ) ;        // target
    EXPECT_EQ( cnot.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CNOT_check( 1 , 0 , 0 , 0 ,
                                                 0 , 1 , 0 , 0 ,
                                                 0 , 0 , 0 , 1 ,
                                                 0 , 0 , 1 , 0 ) ;
    EXPECT_TRUE( cnot.matrix() == CNOT_check ) ;

    // qubit
    EXPECT_EQ( cnot.qubit() , 0 ) ;

    // qubits
    auto qubits = cnot.qubits() ;
    EXPECT_EQ( qubits.size() , 2 ) ;
    EXPECT_EQ( qubits[0] , 0 ) ;
    EXPECT_EQ( qubits[1] , 1 ) ;
    int qnew[] = { 5 , 3 } ;
    cnot.setQubits( &qnew[0] ) ;
    EXPECT_EQ( cnot.qubits()[0] , 3 ) ;
    EXPECT_EQ( cnot.qubits()[1] , 5 ) ;
    qnew[0] = 0 ;
    qnew[1] = 1 ;
    cnot.setQubits( &qnew[0] ) ;

    // print
    cnot.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cnot.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cx q[0], q[1];\n" ) ;
    std::cout << qasm.str() ;

    // gate
    qclab::qgates::PauliX< T >  X ;
    EXPECT_TRUE( *cnot.gate() == X ) ;
    EXPECT_TRUE( cnot.gate()->matrix() == X.matrix() ) ;

    // operators == and !=
    EXPECT_TRUE(  cnot != X ) ;
    EXPECT_FALSE( cnot == X ) ;
    qclab::qgates::CNOT< T >  cnot2 ;
    EXPECT_TRUE(  cnot == cnot2 ) ;
    EXPECT_FALSE( cnot != cnot2 ) ;

    // setControl, setTarget, setControlState
    cnot.setControl( 3 ) ;
    EXPECT_EQ( cnot.control() , 3 ) ;
    cnot.setTarget( 5 ) ;
    EXPECT_EQ( cnot.target() , 5 ) ;
    EXPECT_TRUE(  cnot == cnot2 ) ;
    EXPECT_FALSE( cnot != cnot2 ) ;

    cnot.setControl( 4 ) ;
    EXPECT_EQ( cnot.control() , 4 ) ;
    cnot.setTarget( 1 ) ;
    EXPECT_EQ( cnot.target() , 1 ) ;
    EXPECT_TRUE(  cnot != cnot2 ) ;
    EXPECT_FALSE( cnot == cnot2 ) ;

    qubits[0] = 1 ;
    qubits[1] = 2 ;
    cnot.setQubits( &qubits[0] ) ;
    EXPECT_EQ( cnot.control() , 1 ) ;
    EXPECT_EQ( cnot.target() , 2 ) ;
    EXPECT_TRUE(  cnot == cnot2 ) ;
    EXPECT_FALSE( cnot != cnot2 ) ;

    cnot.setControl( 0 ) ;
    EXPECT_EQ( cnot.control() , 0 ) ;
    cnot.setTarget( 1 ) ;
    EXPECT_EQ( cnot.target() , 1 ) ;
    cnot.setControlState( 0 ) ;
    EXPECT_EQ( cnot.controlState() , 0 ) ;
    EXPECT_TRUE(  cnot != cnot2 ) ;
    EXPECT_FALSE( cnot == cnot2 ) ;
  }

  {
    qclab::qgates::CNOT< T >  cnot01( 0 , 1 ) ;

    EXPECT_EQ( cnot01.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cnot01.fixed() ) ;           // fixed
    EXPECT_TRUE( cnot01.controlled() ) ;      // controlled
    EXPECT_EQ( cnot01.control() , 0 ) ;       // control
    EXPECT_EQ( cnot01.target() , 1 ) ;        // target
    EXPECT_EQ( cnot01.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CNOT_check( 1 , 0 , 0 , 0 ,
                                                 0 , 1 , 0 , 0 ,
                                                 0 , 0 , 0 , 1 ,
                                                 0 , 0 , 1 , 0 ) ;
    EXPECT_TRUE( cnot01.matrix() == CNOT_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cnot01.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cx q[0], q[1];\n" ) ;
  }

  {
    qclab::qgates::CNOT< T >  cnot35( 3 , 5 ) ;

    EXPECT_EQ( cnot35.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cnot35.fixed() ) ;           // fixed
    EXPECT_TRUE( cnot35.controlled() ) ;      // controlled
    EXPECT_EQ( cnot35.control() , 3 ) ;       // control
    EXPECT_EQ( cnot35.target() , 5 ) ;        // target
    EXPECT_EQ( cnot35.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CNOT_check( 1 , 0 , 0 , 0 ,
                                                 0 , 1 , 0 , 0 ,
                                                 0 , 0 , 0 , 1 ,
                                                 0 , 0 , 1 , 0 ) ;
    EXPECT_TRUE( cnot35.matrix() == CNOT_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cnot35.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cx q[3], q[5];\n" ) ;
  }

  {
    qclab::qgates::CNOT< T >  cnot10( 1 , 0 ) ;

    EXPECT_EQ( cnot10.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cnot10.fixed() ) ;           // fixed
    EXPECT_TRUE( cnot10.controlled() ) ;      // controlled
    EXPECT_EQ( cnot10.control() , 1 ) ;       // control
    EXPECT_EQ( cnot10.target() , 0 ) ;        // target
    EXPECT_EQ( cnot10.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CNOT_check( 1 , 0 , 0 , 0 ,
                                                 0 , 0 , 0 , 1 ,
                                                 0 , 0 , 1 , 0 ,
                                                 0 , 1 , 0 , 0 ) ;
    EXPECT_TRUE( cnot10.matrix() == CNOT_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cnot10.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cx q[1], q[0];\n" ) ;
  }

  {
    qclab::qgates::CNOT< T >  cnot53( 5 , 3 ) ;

    EXPECT_EQ( cnot53.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cnot53.fixed() ) ;           // fixed
    EXPECT_TRUE( cnot53.controlled() ) ;      // controlled
    EXPECT_EQ( cnot53.control() , 5 ) ;       // control
    EXPECT_EQ( cnot53.target() , 3 ) ;        // target
    EXPECT_EQ( cnot53.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CNOT_check( 1 , 0 , 0 , 0 ,
                                                 0 , 0 , 0 , 1 ,
                                                 0 , 0 , 1 , 0 ,
                                                 0 , 1 , 0 , 0 ) ;
    EXPECT_TRUE( cnot53.matrix() == CNOT_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cnot53.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cx q[5], q[3];\n" ) ;
  }


  //
  // CNOTs |0> controlled
  //
  {
    qclab::qgates::CNOT< T >  cnot01( 0 , 1 , 0 ) ;

    EXPECT_EQ( cnot01.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cnot01.fixed() ) ;           // fixed
    EXPECT_TRUE( cnot01.controlled() ) ;      // controlled
    EXPECT_EQ( cnot01.control() , 0 ) ;       // control
    EXPECT_EQ( cnot01.target() , 1 ) ;        // target
    EXPECT_EQ( cnot01.controlState() , 0 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CNOT_check( 0 , 1 , 0 , 0 ,
                                                 1 , 0 , 0 , 0 ,
                                                 0 , 0 , 1 , 0 ,
                                                 0 , 0 , 0 , 1 ) ;
    EXPECT_TRUE( cnot01.matrix() == CNOT_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cnot01.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[0];\ncx q[0], q[1];\nx q[0];\n" ) ;
  }

  {
    qclab::qgates::CNOT< T >  cnot10( 1 , 0 , 0 ) ;

    EXPECT_EQ( cnot10.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cnot10.fixed() ) ;           // fixed
    EXPECT_TRUE( cnot10.controlled() ) ;      // controlled
    EXPECT_EQ( cnot10.control() , 1 ) ;       // control
    EXPECT_EQ( cnot10.target() , 0 ) ;        // target
    EXPECT_EQ( cnot10.controlState() , 0 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CNOT_check( 0 , 0 , 1 , 0 ,
                                                 0 , 1 , 0 , 0 ,
                                                 1 , 0 , 0 , 0 ,
                                                 0 , 0 , 0 , 1 ) ;
    EXPECT_TRUE( cnot10.matrix() == CNOT_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cnot10.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[1];\ncx q[1], q[0];\nx q[1];\n" ) ;
  }

}


/*
 * float
 */
TEST( qclab_qgates_CNOT , float ) {
  test_qclab_qgates_CNOT< float >() ;
}

/*
 * double
 */
TEST( qclab_qgates_CNOT , double ) {
  test_qclab_qgates_CNOT< double >() ;
}


/*
 * complex float
 */
TEST( qclab_qgates_CNOT , complex_float ) {
  test_qclab_qgates_CNOT< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_CNOT , complex_double ) {
  test_qclab_qgates_CNOT< std::complex< double > >() ;
}

