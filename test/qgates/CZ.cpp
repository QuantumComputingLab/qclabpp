#include <gtest/gtest.h>
#include "qclab/qgates/CZ.hpp"

template <typename T>
void test_qclab_qgates_CZ() {

  //
  // CZs |1> controlled
  //
  {
    qclab::qgates::CZ< T >  cz ;

    EXPECT_EQ( cz.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cz.fixed() ) ;           // fixed
    EXPECT_TRUE( cz.controlled() ) ;      // controlled
    EXPECT_EQ( cz.control() , 0 ) ;       // control
    EXPECT_EQ( cz.target() , 1 ) ;        // target
    EXPECT_EQ( cz.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CZ_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 0 , 0 ,-1 ) ;
    EXPECT_TRUE( cz.matrix() == CZ_check ) ;

    // qubit
    EXPECT_EQ( cz.qubit() , 0 ) ;

    // qubits
    auto qubits = cz.qubits() ;
    EXPECT_EQ( qubits.size() , 2 ) ;
    EXPECT_EQ( qubits[0] , 0 ) ;
    EXPECT_EQ( qubits[1] , 1 ) ;
    int qnew[] = { 5 , 3 } ;
    cz.setQubits( &qnew[0] ) ;
    EXPECT_EQ( cz.qubits()[0] , 3 ) ;
    EXPECT_EQ( cz.qubits()[1] , 5 ) ;
    qnew[0] = 0 ;
    qnew[1] = 1 ;
    cz.setQubits( &qnew[0] ) ;

    // print
    cz.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cz.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cz q[0], q[1];\n" ) ;
    std::cout << qasm.str() ;

    // gate
    qclab::qgates::PauliZ< T >  Z ;
    EXPECT_TRUE( *cz.gate() == Z ) ;
    EXPECT_TRUE( cz.gate()->matrix() == Z.matrix() ) ;

    // operators == and !=
    EXPECT_TRUE(  cz != Z ) ;
    EXPECT_FALSE( cz == Z ) ;
    qclab::qgates::CZ< T >  cz2 ;
    EXPECT_TRUE(  cz == cz2 ) ;
    EXPECT_FALSE( cz != cz2 ) ;

    // setControl, setTarget, setControlState
    cz.setControl( 3 ) ;
    EXPECT_EQ( cz.control() , 3 ) ;
    cz.setTarget( 5 ) ;
    EXPECT_EQ( cz.target() , 5 ) ;
    EXPECT_TRUE(  cz == cz2 ) ;
    EXPECT_FALSE( cz != cz2 ) ;

    cz.setControl( 4 ) ;
    EXPECT_EQ( cz.control() , 4 ) ;
    cz.setTarget( 1 ) ;
    EXPECT_EQ( cz.target() , 1 ) ;
    EXPECT_TRUE(  cz != cz2 ) ;
    EXPECT_FALSE( cz == cz2 ) ;

    qubits[0] = 1 ;
    qubits[1] = 2 ;
    cz.setQubits( &qubits[0] ) ;
    EXPECT_EQ( cz.control() , 1 ) ;
    EXPECT_EQ( cz.target() , 2 ) ;
    EXPECT_TRUE(  cz == cz2 ) ;
    EXPECT_FALSE( cz != cz2 ) ;

    cz.setControl( 0 ) ;
    EXPECT_EQ( cz.control() , 0 ) ;
    cz.setTarget( 1 ) ;
    EXPECT_EQ( cz.target() , 1 ) ;
    cz.setControlState( 0 ) ;
    EXPECT_EQ( cz.controlState() , 0 ) ;
    EXPECT_TRUE(  cz != cz2 ) ;
    EXPECT_FALSE( cz == cz2 ) ;
  }

  {
    qclab::qgates::CZ< T >  cz01( 0 , 1 ) ;

    EXPECT_EQ( cz01.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cz01.fixed() ) ;           // fixed
    EXPECT_TRUE( cz01.controlled() ) ;      // controlled
    EXPECT_EQ( cz01.control() , 0 ) ;       // control
    EXPECT_EQ( cz01.target() , 1 ) ;        // target
    EXPECT_EQ( cz01.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CZ_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 0 , 0 ,-1 ) ;
    EXPECT_TRUE( cz01.matrix() == CZ_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cz01.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cz q[0], q[1];\n" ) ;
  }

  {
    qclab::qgates::CZ< T >  cz35( 3 , 5 ) ;

    EXPECT_EQ( cz35.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cz35.fixed() ) ;           // fixed
    EXPECT_TRUE( cz35.controlled() ) ;      // controlled
    EXPECT_EQ( cz35.control() , 3 ) ;       // control
    EXPECT_EQ( cz35.target() , 5 ) ;        // target
    EXPECT_EQ( cz35.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CZ_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 0 , 0 ,-1 ) ;
    EXPECT_TRUE( cz35.matrix() == CZ_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cz35.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cz q[3], q[5];\n" ) ;
  }

  {
    qclab::qgates::CZ< T >  cz10( 1 , 0 ) ;

    EXPECT_EQ( cz10.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cz10.fixed() ) ;           // fixed
    EXPECT_TRUE( cz10.controlled() ) ;      // controlled
    EXPECT_EQ( cz10.control() , 1 ) ;       // control
    EXPECT_EQ( cz10.target() , 0 ) ;        // target
    EXPECT_EQ( cz10.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CZ_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 0 , 0 ,-1 ) ;
    EXPECT_TRUE( cz10.matrix() == CZ_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cz10.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cz q[1], q[0];\n" ) ;
  }

  {
    qclab::qgates::CZ< T >  cz53( 5 , 3 ) ;

    EXPECT_EQ( cz53.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cz53.fixed() ) ;           // fixed
    EXPECT_TRUE( cz53.controlled() ) ;      // controlled
    EXPECT_EQ( cz53.control() , 5 ) ;       // control
    EXPECT_EQ( cz53.target() , 3 ) ;        // target
    EXPECT_EQ( cz53.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CZ_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 0 , 0 ,-1 ) ;
    EXPECT_TRUE( cz53.matrix() == CZ_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cz53.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cz q[5], q[3];\n" ) ;
  }


  //
  // CZs |0> controlled
  //
  {
    qclab::qgates::CZ< T >  cz01( 0 , 1 , 0 ) ;

    EXPECT_EQ( cz01.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cz01.fixed() ) ;           // fixed
    EXPECT_TRUE( cz01.controlled() ) ;      // controlled
    EXPECT_EQ( cz01.control() , 0 ) ;       // control
    EXPECT_EQ( cz01.target() , 1 ) ;        // target
    EXPECT_EQ( cz01.controlState() , 0 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CZ_check( 1 , 0 , 0 , 0 ,
                                               0 ,-1 , 0 , 0 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 0 , 0 , 1 ) ;
    EXPECT_TRUE( cz01.matrix() == CZ_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cz01.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[0];\ncz q[0], q[1];\nx q[0];\n" ) ;
  }

  {
    qclab::qgates::CZ< T >  cz10( 1 , 0 , 0 ) ;

    EXPECT_EQ( cz10.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cz10.fixed() ) ;           // fixed
    EXPECT_TRUE( cz10.controlled() ) ;      // controlled
    EXPECT_EQ( cz10.control() , 1 ) ;       // control
    EXPECT_EQ( cz10.target() , 0 ) ;        // target
    EXPECT_EQ( cz10.controlState() , 0 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CZ_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 ,-1 , 0 ,
                                               0 , 0 , 0 , 1 ) ;
    EXPECT_TRUE( cz10.matrix() == CZ_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cz10.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[1];\ncz q[1], q[0];\nx q[1];\n" ) ;
  }

}


/*
 * float
 */
TEST( qclab_qgates_CZ , float ) {
  test_qclab_qgates_CZ< float >() ;
}

/*
 * double
 */
TEST( qclab_qgates_CZ , double ) {
  test_qclab_qgates_CZ< double >() ;
}


/*
 * complex float
 */
TEST( qclab_qgates_CZ , complex_float ) {
  test_qclab_qgates_CZ< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_CZ , complex_double ) {
  test_qclab_qgates_CZ< std::complex< double > >() ;
}

