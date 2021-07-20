#include <gtest/gtest.h>
#include "qclab/qgates/CX.hpp"

template <typename T>
void test_qclab_qgates_CX() {

  //
  // CXs |1> controlled
  //
  {
    qclab::qgates::CX< T >  cx ;

    EXPECT_EQ( cx.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cx.fixed() ) ;           // fixed
    EXPECT_TRUE( cx.controlled() ) ;      // controlled
    EXPECT_EQ( cx.control() , 0 ) ;       // control
    EXPECT_EQ( cx.target() , 1 ) ;        // target
    EXPECT_EQ( cx.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CX_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 0 , 1 ,
                                               0 , 0 , 1 , 0 ) ;
    EXPECT_TRUE( cx.matrix() == CX_check ) ;

    // qubit
    EXPECT_EQ( cx.qubit() , 0 ) ;

    // qubits
    auto qubits = cx.qubits() ;
    EXPECT_EQ( qubits.size() , 2 ) ;
    EXPECT_EQ( qubits[0] , 0 ) ;
    EXPECT_EQ( qubits[1] , 1 ) ;
    int qnew[] = { 5 , 3 } ;
    cx.setQubits( &qnew[0] ) ;
    EXPECT_EQ( cx.qubits()[0] , 3 ) ;
    EXPECT_EQ( cx.qubits()[1] , 5 ) ;
    qnew[0] = 0 ;
    qnew[1] = 1 ;
    cx.setQubits( &qnew[0] ) ;

    // print
    cx.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cx.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cx q[0], q[1];\n" ) ;
    std::cout << qasm.str() ;

    // gate
    qclab::qgates::PauliX< T >  X ;
    EXPECT_TRUE( *cx.gate() == X ) ;
    EXPECT_TRUE( cx.gate()->matrix() == X.matrix() ) ;

    // operators == and !=
    EXPECT_TRUE(  cx != X ) ;
    EXPECT_FALSE( cx == X ) ;
    qclab::qgates::CX< T >  cx2 ;
    EXPECT_TRUE(  cx == cx2 ) ;
    EXPECT_FALSE( cx != cx2 ) ;

    // setControl, setTarget, setControlState
    cx.setControl( 3 ) ;
    EXPECT_EQ( cx.control() , 3 ) ;
    cx.setTarget( 5 ) ;
    EXPECT_EQ( cx.target() , 5 ) ;
    EXPECT_TRUE(  cx == cx2 ) ;
    EXPECT_FALSE( cx != cx2 ) ;

    cx.setControl( 4 ) ;
    EXPECT_EQ( cx.control() , 4 ) ;
    cx.setTarget( 1 ) ;
    EXPECT_EQ( cx.target() , 1 ) ;
    EXPECT_TRUE(  cx != cx2 ) ;
    EXPECT_FALSE( cx == cx2 ) ;

    qubits[0] = 1 ;
    qubits[1] = 2 ;
    cx.setQubits( &qubits[0] ) ;
    EXPECT_EQ( cx.control() , 1 ) ;
    EXPECT_EQ( cx.target() , 2 ) ;
    EXPECT_TRUE(  cx == cx2 ) ;
    EXPECT_FALSE( cx != cx2 ) ;

    cx.setControl( 0 ) ;
    EXPECT_EQ( cx.control() , 0 ) ;
    cx.setTarget( 1 ) ;
    EXPECT_EQ( cx.target() , 1 ) ;
    cx.setControlState( 0 ) ;
    EXPECT_EQ( cx.controlState() , 0 ) ;
    EXPECT_TRUE(  cx != cx2 ) ;
    EXPECT_FALSE( cx == cx2 ) ;
  }

  {
    qclab::qgates::CX< T >  cx01( 0 , 1 ) ;

    EXPECT_EQ( cx01.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cx01.fixed() ) ;           // fixed
    EXPECT_TRUE( cx01.controlled() ) ;      // controlled
    EXPECT_EQ( cx01.control() , 0 ) ;       // control
    EXPECT_EQ( cx01.target() , 1 ) ;        // target
    EXPECT_EQ( cx01.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CX_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 0 , 1 ,
                                               0 , 0 , 1 , 0 ) ;
    EXPECT_TRUE( cx01.matrix() == CX_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cx01.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cx q[0], q[1];\n" ) ;
  }

  {
    qclab::qgates::CX< T >  cx35( 3 , 5 ) ;

    EXPECT_EQ( cx35.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cx35.fixed() ) ;           // fixed
    EXPECT_TRUE( cx35.controlled() ) ;      // controlled
    EXPECT_EQ( cx35.control() , 3 ) ;       // control
    EXPECT_EQ( cx35.target() , 5 ) ;        // target
    EXPECT_EQ( cx35.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CX_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 0 , 1 ,
                                               0 , 0 , 1 , 0 ) ;
    EXPECT_TRUE( cx35.matrix() == CX_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cx35.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cx q[3], q[5];\n" ) ;
  }

  {
    qclab::qgates::CX< T >  cx10( 1 , 0 ) ;

    EXPECT_EQ( cx10.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cx10.fixed() ) ;           // fixed
    EXPECT_TRUE( cx10.controlled() ) ;      // controlled
    EXPECT_EQ( cx10.control() , 1 ) ;       // control
    EXPECT_EQ( cx10.target() , 0 ) ;        // target
    EXPECT_EQ( cx10.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CX_check( 1 , 0 , 0 , 0 ,
                                               0 , 0 , 0 , 1 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 1 , 0 , 0 ) ;
    EXPECT_TRUE( cx10.matrix() == CX_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cx10.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cx q[1], q[0];\n" ) ;
  }

  {
    qclab::qgates::CX< T >  cx53( 5 , 3 ) ;

    EXPECT_EQ( cx53.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cx53.fixed() ) ;           // fixed
    EXPECT_TRUE( cx53.controlled() ) ;      // controlled
    EXPECT_EQ( cx53.control() , 5 ) ;       // control
    EXPECT_EQ( cx53.target() , 3 ) ;        // target
    EXPECT_EQ( cx53.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CX_check( 1 , 0 , 0 , 0 ,
                                               0 , 0 , 0 , 1 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 1 , 0 , 0 ) ;
    EXPECT_TRUE( cx53.matrix() == CX_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cx53.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cx q[5], q[3];\n" ) ;
  }


  //
  // CXs |0> controlled
  //
  {
    qclab::qgates::CX< T >  cx01( 0 , 1 , 0 ) ;

    EXPECT_EQ( cx01.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cx01.fixed() ) ;           // fixed
    EXPECT_TRUE( cx01.controlled() ) ;      // controlled
    EXPECT_EQ( cx01.control() , 0 ) ;       // control
    EXPECT_EQ( cx01.target() , 1 ) ;        // target
    EXPECT_EQ( cx01.controlState() , 0 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CX_check( 0 , 1 , 0 , 0 ,
                                               1 , 0 , 0 , 0 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 0 , 0 , 1 ) ;
    EXPECT_TRUE( cx01.matrix() == CX_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cx01.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[0];\ncx q[0], q[1];\nx q[0];\n" ) ;
  }

  {
    qclab::qgates::CX< T >  cx10( 1 , 0 , 0 ) ;

    EXPECT_EQ( cx10.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cx10.fixed() ) ;           // fixed
    EXPECT_TRUE( cx10.controlled() ) ;      // controlled
    EXPECT_EQ( cx10.control() , 1 ) ;       // control
    EXPECT_EQ( cx10.target() , 0 ) ;        // target
    EXPECT_EQ( cx10.controlState() , 0 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CX_check( 0 , 0 , 1 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               1 , 0 , 0 , 0 ,
                                               0 , 0 , 0 , 1 ) ;
    EXPECT_TRUE( cx10.matrix() == CX_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cx10.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[1];\ncx q[1], q[0];\nx q[1];\n" ) ;
  }

}


/*
 * float
 */
TEST( qclab_qgates_CX , float ) {
  test_qclab_qgates_CX< float >() ;
}

/*
 * double
 */
TEST( qclab_qgates_CX , double ) {
  test_qclab_qgates_CX< double >() ;
}


/*
 * complex float
 */
TEST( qclab_qgates_CX , complex_float ) {
  test_qclab_qgates_CX< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_CX , complex_double ) {
  test_qclab_qgates_CX< std::complex< double > >() ;
}

