#include <gtest/gtest.h>
#include "qclab/qgates/SWAP.hpp"

template <typename T>
void test_qclab_qgates_SWAP() {

  const auto I1 = qclab::dense::eye< T >(  2 ) ;
  const auto I2 = qclab::dense::eye< T >(  4 ) ;
  const auto I3 = qclab::dense::eye< T >(  8 ) ;
  const auto I4 = qclab::dense::eye< T >( 16 ) ;

  {
    qclab::qgates::SWAP< T >  swap ;

    EXPECT_EQ( swap.nbQubits() , 2 ) ;   // nbQubits
    EXPECT_TRUE( swap.fixed() ) ;        // fixed
    EXPECT_FALSE( swap.controlled() ) ;  // controlled

    // qubit
    EXPECT_EQ( swap.qubit() , 0 ) ;

    // qubits
    auto qubits = swap.qubits() ;
    EXPECT_EQ( qubits.size() , 2 ) ;
    EXPECT_EQ( qubits[0] , 0 ) ;
    EXPECT_EQ( qubits[1] , 1 ) ;
    int qnew[] = { 3 , 5 } ;
    swap.setQubits( &qnew[0] ) ;
    EXPECT_EQ( swap.qubits()[0] , 3 ) ;
    EXPECT_EQ( swap.qubits()[1] , 5 ) ;
    qnew[0] = 0 ;
    qnew[1] = 1 ;
    swap.setQubits( &qnew[0] ) ;

    // matrix
    qclab::dense::SquareMatrix< T >  SWAP_check( 1 , 0 , 0 , 0 ,
                                                 0 , 0 , 1 , 0 ,
                                                 0 , 1 , 0 , 0 ,
                                                 0 , 0 , 0 , 1 ) ;
    EXPECT_TRUE( swap.matrix() == SWAP_check ) ;

    // print
    swap.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( swap.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "swap q[0], q[1];\n" ) ;
    std::cout << qasm.str() ;

    // operators == and !=
    qclab::qgates::SWAP< T >  swap2( 2 , 4 ) ;
    EXPECT_TRUE(  swap == swap2 ) ;
    EXPECT_FALSE( swap != swap2 ) ;
  }

  {
    qclab::qgates::SWAP< T >  swap( 3 , 5 ) ;

    EXPECT_EQ( swap.nbQubits() , 2 ) ;   // nbQubits
    EXPECT_TRUE( swap.fixed() ) ;        // fixed
    EXPECT_FALSE( swap.controlled() ) ;  // controlled
    EXPECT_EQ( swap.qubits()[0] , 3 ) ;  // qubit0
    EXPECT_EQ( swap.qubits()[1] , 5 ) ;  // qubit1
  }

  {
    int qubits[] = { 3 , 5 } ;
    qclab::qgates::SWAP< T >  swap( &qubits[0] ) ;

    EXPECT_EQ( swap.nbQubits() , 2 ) ;   // nbQubits
    EXPECT_TRUE( swap.fixed() ) ;        // fixed
    EXPECT_FALSE( swap.controlled() ) ;  // controlled
    EXPECT_EQ( swap.qubits()[0] , 3 ) ;  // qubit0
    EXPECT_EQ( swap.qubits()[1] , 5 ) ;  // qubit1
  }

  {
    qclab::qgates::SWAP< T >  swap( 0 , 1 ) ;

    // apply (2 qubits)
    auto mat2 = I2 ;
    swap.apply( qclab::Side::Left , qclab::Op::NoTrans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == swap.matrix() ) ;

    qclab::dense::SquareMatrix< T >  mat(  1 ,  2 ,  3 ,  4 ,
                                           5 ,  6 ,  7 ,  8 ,
                                           9 , 10 , 11 , 12 ,
                                          13 , 14 , 15 , 16 ) ;

    mat2 = mat ;
    swap.apply( qclab::Side::Left , qclab::Op::NoTrans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == mat * swap.matrix() ) ;

    mat2 = mat ;
    swap.apply( qclab::Side::Right , qclab::Op::NoTrans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == swap.matrix() * mat ) ;

    // apply (3 qubits)
    auto mat3 = I3 ;
    swap.apply( qclab::Side::Left , qclab::Op::NoTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( swap.matrix() , I1 ) ) ;

    int qnew[] = { 1 , 2 } ;
    swap.setQubits( &qnew[0] ) ;
    mat3 = I3 ;
    swap.apply( qclab::Side::Left , qclab::Op::NoTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( I1 , swap.matrix() ) ) ;
  }

  {
    qclab::qgates::SWAP< T >  swap( 0 , 2 ) ;

    auto check = qclab::dense::zeros< T >( 8 ) ;
    check(0,0) = 1 ;
    check(4,1) = 1 ;
    check(2,2) = 1 ;
    check(6,3) = 1 ;
    check(1,4) = 1 ;
    check(5,5) = 1 ;
    check(3,6) = 1 ;
    check(7,7) = 1 ;

    // apply (3 qubits)
    auto mat3 = I3 ;
    swap.apply( qclab::Side::Left , qclab::Op::NoTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == check ) ;
  }

}


/*
 * float
 */
TEST( qclab_qgates_SWAP , float ) {
  test_qclab_qgates_SWAP< float >() ;
}

/*
 * double
 */
TEST( qclab_qgates_SWAP , double ) {
  test_qclab_qgates_SWAP< double >() ;
}


/*
 * complex float
 */
TEST( qclab_qgates_SWAP , complex_float ) {
  test_qclab_qgates_SWAP< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_SWAP , complex_double ) {
  test_qclab_qgates_SWAP< std::complex< double > >() ;
}

