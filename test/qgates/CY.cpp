#include <gtest/gtest.h>
#include "qclab/qgates/CY.hpp"

template <typename T>
void test_qclab_qgates_CY() {

  //
  // CYs |1> controlled
  //
  {
    qclab::qgates::CY< T >  cy ;

    EXPECT_EQ( cy.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cy.fixed() ) ;           // fixed
    EXPECT_TRUE( cy.controlled() ) ;      // controlled
    EXPECT_EQ( cy.control() , 0 ) ;       // control
    EXPECT_EQ( cy.target() , 1 ) ;        // target
    EXPECT_EQ( cy.controlState() , 1 ) ;  // controlState

    // matrix
    const T i(0,1) ;
    qclab::dense::SquareMatrix< T >  CY_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 0 ,-i ,
                                               0 , 0 , i , 0 ) ;
    EXPECT_TRUE( cy.matrix() == CY_check ) ;

    // qubit
    EXPECT_EQ( cy.qubit() , 0 ) ;

    // qubits
    auto qubits = cy.qubits() ;
    EXPECT_EQ( qubits.size() , 2 ) ;
    EXPECT_EQ( qubits[0] , 0 ) ;
    EXPECT_EQ( qubits[1] , 1 ) ;
    int qnew[] = { 5 , 3 } ;
    cy.setQubits( &qnew[0] ) ;
    EXPECT_EQ( cy.qubits()[0] , 3 ) ;
    EXPECT_EQ( cy.qubits()[1] , 5 ) ;
    qnew[0] = 0 ;
    qnew[1] = 1 ;
    cy.setQubits( &qnew[0] ) ;

    // print
    cy.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cy.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cy q[0], q[1];\n" ) ;
    std::cout << qasm.str() ;

    // gate
    qclab::qgates::PauliY< T >  Y ;
    EXPECT_TRUE( *cy.gate() == Y ) ;
    EXPECT_TRUE( cy.gate()->matrix() == Y.matrix() ) ;

    // operators == and !=
    EXPECT_TRUE(  cy != Y ) ;
    EXPECT_FALSE( cy == Y ) ;
    qclab::qgates::CY< T >  cy2 ;
    EXPECT_TRUE(  cy == cy2 ) ;
    EXPECT_FALSE( cy != cy2 ) ;

    // setControl, setTarget, setControlState
    cy.setControl( 3 ) ;
    EXPECT_EQ( cy.control() , 3 ) ;
    cy.setTarget( 5 ) ;
    EXPECT_EQ( cy.target() , 5 ) ;
    EXPECT_TRUE(  cy == cy2 ) ;
    EXPECT_FALSE( cy != cy2 ) ;

    cy.setControl( 4 ) ;
    EXPECT_EQ( cy.control() , 4 ) ;
    cy.setTarget( 1 ) ;
    EXPECT_EQ( cy.target() , 1 ) ;
    EXPECT_TRUE(  cy != cy2 ) ;
    EXPECT_FALSE( cy == cy2 ) ;

    qubits[0] = 1 ;
    qubits[1] = 2 ;
    cy.setQubits( &qubits[0] ) ;
    EXPECT_EQ( cy.control() , 1 ) ;
    EXPECT_EQ( cy.target() , 2 ) ;
    EXPECT_TRUE(  cy == cy2 ) ;
    EXPECT_FALSE( cy != cy2 ) ;

    cy.setControl( 0 ) ;
    EXPECT_EQ( cy.control() , 0 ) ;
    cy.setTarget( 1 ) ;
    EXPECT_EQ( cy.target() , 1 ) ;
    cy.setControlState( 0 ) ;
    EXPECT_EQ( cy.controlState() , 0 ) ;
    EXPECT_TRUE(  cy != cy2 ) ;
    EXPECT_FALSE( cy == cy2 ) ;
  }

  {
    qclab::qgates::CY< T >  cy01( 0 , 1 ) ;

    EXPECT_EQ( cy01.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cy01.fixed() ) ;           // fixed
    EXPECT_TRUE( cy01.controlled() ) ;      // controlled
    EXPECT_EQ( cy01.control() , 0 ) ;       // control
    EXPECT_EQ( cy01.target() , 1 ) ;        // target
    EXPECT_EQ( cy01.controlState() , 1 ) ;  // controlState

    // matrix
    const T i(0,1) ;
    qclab::dense::SquareMatrix< T >  CY_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 0 ,-i ,
                                               0 , 0 , i , 0 ) ;
    EXPECT_TRUE( cy01.matrix() == CY_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cy01.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cy q[0], q[1];\n" ) ;
  }

  {
    qclab::qgates::CY< T >  cy35( 3 , 5 ) ;

    EXPECT_EQ( cy35.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cy35.fixed() ) ;           // fixed
    EXPECT_TRUE( cy35.controlled() ) ;      // controlled
    EXPECT_EQ( cy35.control() , 3 ) ;       // control
    EXPECT_EQ( cy35.target() , 5 ) ;        // target
    EXPECT_EQ( cy35.controlState() , 1 ) ;  // controlState

    // matrix
    const T i(0,1) ;
    qclab::dense::SquareMatrix< T >  CY_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 0 ,-i ,
                                               0 , 0 , i , 0 ) ;
    EXPECT_TRUE( cy35.matrix() == CY_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cy35.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cy q[3], q[5];\n" ) ;
  }

  {
    qclab::qgates::CY< T >  cy10( 1 , 0 ) ;

    EXPECT_EQ( cy10.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cy10.fixed() ) ;           // fixed
    EXPECT_TRUE( cy10.controlled() ) ;      // controlled
    EXPECT_EQ( cy10.control() , 1 ) ;       // control
    EXPECT_EQ( cy10.target() , 0 ) ;        // target
    EXPECT_EQ( cy10.controlState() , 1 ) ;  // controlState

    // matrix
    const T i(0,1) ;
    qclab::dense::SquareMatrix< T >  CY_check( 1 , 0 , 0 , 0 ,
                                               0 , 0 , 0 ,-i ,
                                               0 , 0 , 1 , 0 ,
                                               0 , i , 0 , 0 ) ;
    EXPECT_TRUE( cy10.matrix() == CY_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cy10.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cy q[1], q[0];\n" ) ;
  }

  {
    qclab::qgates::CY< T >  cy53( 5 , 3 ) ;

    EXPECT_EQ( cy53.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cy53.fixed() ) ;           // fixed
    EXPECT_TRUE( cy53.controlled() ) ;      // controlled
    EXPECT_EQ( cy53.control() , 5 ) ;       // control
    EXPECT_EQ( cy53.target() , 3 ) ;        // target
    EXPECT_EQ( cy53.controlState() , 1 ) ;  // controlState

    // matrix
    const T i(0,1) ;
    qclab::dense::SquareMatrix< T >  CY_check( 1 , 0 , 0 , 0 ,
                                               0 , 0 , 0 ,-i ,
                                               0 , 0 , 1 , 0 ,
                                               0 , i , 0 , 0 ) ;
    EXPECT_TRUE( cy53.matrix() == CY_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cy53.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cy q[5], q[3];\n" ) ;
  }


  //
  // CYs |0> controlled
  //
  {
    qclab::qgates::CY< T >  cy01( 0 , 1 , 0 ) ;

    EXPECT_EQ( cy01.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cy01.fixed() ) ;           // fixed
    EXPECT_TRUE( cy01.controlled() ) ;      // controlled
    EXPECT_EQ( cy01.control() , 0 ) ;       // control
    EXPECT_EQ( cy01.target() , 1 ) ;        // target
    EXPECT_EQ( cy01.controlState() , 0 ) ;  // controlState

    // matrix
    const T i(0,1) ;
    qclab::dense::SquareMatrix< T >  CY_check( 0 ,-i , 0 , 0 ,
                                               i , 0 , 0 , 0 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 0 , 0 , 1 ) ;
    EXPECT_TRUE( cy01.matrix() == CY_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cy01.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[0];\ncy q[0], q[1];\nx q[0];\n" ) ;
  }

  {
    qclab::qgates::CY< T >  cy10( 1 , 0 , 0 ) ;

    EXPECT_EQ( cy10.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cy10.fixed() ) ;           // fixed
    EXPECT_TRUE( cy10.controlled() ) ;      // controlled
    EXPECT_EQ( cy10.control() , 1 ) ;       // control
    EXPECT_EQ( cy10.target() , 0 ) ;        // target
    EXPECT_EQ( cy10.controlState() , 0 ) ;  // controlState

    // matrix
    const T i(0,1) ;
    qclab::dense::SquareMatrix< T >  CY_check( 0 , 0 ,-i , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               i , 0 , 0 , 0 ,
                                               0 , 0 , 0 , 1 ) ;
    EXPECT_TRUE( cy10.matrix() == CY_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cy10.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[1];\ncy q[1], q[0];\nx q[1];\n" ) ;
  }

}


/*
 * complex float
 */
TEST( qclab_qgates_CY , complex_float ) {
  test_qclab_qgates_CY< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_CY , complex_double ) {
  test_qclab_qgates_CY< std::complex< double > >() ;
}

