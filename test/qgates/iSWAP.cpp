#include <gtest/gtest.h>
#include "qgates/iSWAP.hpp"

template <typename M>
void test_qclab_qgates_iSWAP_check( const M m1 , const M m2 ) {

  assert( m1.rows() == m2.rows() ) ;
  assert( m1.cols() == m2.cols() ) ;

  using R = qclab::real_t< typename M::value_type > ;
  const R eps = std::numeric_limits< R >::epsilon() ;
  const R tol = 100 * eps ;

  for ( int j = 0; j < m1.cols(); j++ ) {
    for ( int i = 0; i < m1.rows(); i++ ) {
      EXPECT_NEAR( std::real( m1(i,j) ) , std::real( m2(i,j) ) , tol ) ;
      EXPECT_NEAR( std::imag( m1(i,j) ) , std::imag( m2(i,j) ) , tol ) ;
    }
  }

}


template <typename T>
void test_qclab_qgates_iSWAP() {

  const auto I1 = qclab::dense::eye< T >(  2 ) ;
  const auto I2 = qclab::dense::eye< T >(  4 ) ;
  const auto I3 = qclab::dense::eye< T >(  8 ) ;
  const auto I4 = qclab::dense::eye< T >( 16 ) ;

  {
    qclab::qgates::iSWAP< T >  iswap ;

    EXPECT_EQ( iswap.nbQubits() , 2 ) ;   // nbQubits
    EXPECT_TRUE( iswap.fixed() ) ;        // fixed
    EXPECT_FALSE( iswap.controlled() ) ;  // controlled

    // qubit
    EXPECT_EQ( iswap.qubit() , 0 ) ;

    // qubits
    auto qubits = iswap.qubits() ;
    EXPECT_EQ( qubits.size() , 2 ) ;
    EXPECT_EQ( qubits[0] , 0 ) ;
    EXPECT_EQ( qubits[1] , 1 ) ;
    int qnew[] = { 3 , 5 } ;
    iswap.setQubits( &qnew[0] ) ;
    EXPECT_EQ( iswap.qubits()[0] , 3 ) ;
    EXPECT_EQ( iswap.qubits()[1] , 5 ) ;
    qnew[0] = 0 ;
    qnew[1] = 1 ;
    iswap.setQubits( &qnew[0] ) ;

    // matrix
    const T i(0,1) ;
    qclab::dense::SquareMatrix< T >  iSWAP_check( 1 , 0 , 0 , 0 ,
                                                  0 , 0 , i , 0 ,
                                                  0 , i , 0 , 0 ,
                                                  0 , 0 , 0 , 1 ) ;
    EXPECT_TRUE( iswap.matrix() == iSWAP_check ) ;

    // print
    iswap.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( iswap.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "iswap q[0], q[1];\n" ) ;
    std::cout << qasm.str() ;

    // operators == and !=
    qclab::qgates::iSWAP< T >  iswap2( 2 , 4 ) ;
    EXPECT_TRUE(  iswap == iswap2 ) ;
    EXPECT_FALSE( iswap != iswap2 ) ;
  }

  {
    qclab::qgates::iSWAP< T >  iswap( 3 , 5 ) ;

    EXPECT_EQ( iswap.nbQubits() , 2 ) ;   // nbQubits
    EXPECT_TRUE( iswap.fixed() ) ;        // fixed
    EXPECT_FALSE( iswap.controlled() ) ;  // controlled
    EXPECT_EQ( iswap.qubits()[0] , 3 ) ;  // qubit0
    EXPECT_EQ( iswap.qubits()[1] , 5 ) ;  // qubit1
  }

  {
    int qubits[] = { 3 , 5 } ;
    qclab::qgates::iSWAP< T >  iswap( &qubits[0] ) ;

    EXPECT_EQ( iswap.nbQubits() , 2 ) ;   // nbQubits
    EXPECT_TRUE( iswap.fixed() ) ;        // fixed
    EXPECT_FALSE( iswap.controlled() ) ;  // controlled
    EXPECT_EQ( iswap.qubits()[0] , 3 ) ;  // qubit0
    EXPECT_EQ( iswap.qubits()[1] , 5 ) ;  // qubit1
  }

  {
    qclab::qgates::iSWAP< T >  iswap( 0 , 1 ) ;

    // apply (2 qubits)
    auto mat2 = I2 ;
    iswap.apply( qclab::Side::Left , qclab::Op::NoTrans , 2 , mat2 ) ;
    test_qclab_qgates_iSWAP_check( mat2 , iswap.matrix() ) ;

    qclab::dense::SquareMatrix< T >  mat(  1 ,  2 ,  3 ,  4 ,
                                           5 ,  6 ,  7 ,  8 ,
                                           9 , 10 , 11 , 12 ,
                                          13 , 14 , 15 , 16 ) ;

    mat2 = mat ;
    iswap.apply( qclab::Side::Left , qclab::Op::NoTrans , 2 , mat2 ) ;
    test_qclab_qgates_iSWAP_check( mat2 , mat * iswap.matrix() ) ;

    mat2 = mat ;
    iswap.apply( qclab::Side::Right , qclab::Op::NoTrans , 2 , mat2 ) ;
    test_qclab_qgates_iSWAP_check( mat2 , iswap.matrix() * mat ) ;

    auto iswapH = qclab::dense::conjTrans( iswap.matrix() ) ;

    mat2 = mat ;
    iswap.apply( qclab::Side::Left , qclab::Op::ConjTrans , 2 , mat2 ) ;
    test_qclab_qgates_iSWAP_check( mat2 , mat * iswapH ) ;

    mat2 = mat ;
    iswap.apply( qclab::Side::Right , qclab::Op::ConjTrans , 2 , mat2 ) ;
    test_qclab_qgates_iSWAP_check( mat2 , iswapH * mat ) ;

    // apply (3 qubits)
    auto mat3 = I3 ;
    iswap.apply( qclab::Side::Left , qclab::Op::NoTrans , 3 , mat3 ) ;
    test_qclab_qgates_iSWAP_check( mat3 ,
                                   qclab::dense::kron( iswap.matrix() , I1 ) ) ;

    mat3 = I3 ;
    iswap.apply( qclab::Side::Right , qclab::Op::NoTrans , 3 , mat3 ) ;
    test_qclab_qgates_iSWAP_check( mat3 ,
                                   qclab::dense::kron( iswap.matrix() , I1 ) ) ;

    mat3 = I3 ;
    iswap.apply( qclab::Side::Left , qclab::Op::ConjTrans , 3 , mat3 ) ;
    test_qclab_qgates_iSWAP_check( mat3 , qclab::dense::kron( iswapH , I1 ) ) ;

    mat3 = I3 ;
    iswap.apply( qclab::Side::Right , qclab::Op::ConjTrans , 3 , mat3 ) ;
    test_qclab_qgates_iSWAP_check( mat3 , qclab::dense::kron( iswapH , I1 ) ) ;

    int qnew[] = { 1 , 2 } ;
    iswap.setQubits( &qnew[0] ) ;
    mat3 = I3 ;
    iswap.apply( qclab::Side::Left , qclab::Op::NoTrans , 3 , mat3 ) ;
    test_qclab_qgates_iSWAP_check( mat3 ,
                                   qclab::dense::kron( I1 , iswap.matrix() ) ) ;

    mat3 = I3 ;
    iswap.apply( qclab::Side::Right , qclab::Op::NoTrans , 3 , mat3 ) ;
    test_qclab_qgates_iSWAP_check( mat3 ,
                                   qclab::dense::kron( I1 , iswap.matrix() ) ) ;

    mat3 = I3 ;
    iswap.apply( qclab::Side::Left , qclab::Op::ConjTrans , 3 , mat3 ) ;
    test_qclab_qgates_iSWAP_check( mat3 , qclab::dense::kron( I1 , iswapH ) ) ;

    mat3 = I3 ;
    iswap.apply( qclab::Side::Right , qclab::Op::ConjTrans , 3 , mat3 ) ;
    test_qclab_qgates_iSWAP_check( mat3 , qclab::dense::kron( I1 , iswapH ) ) ;
  }

  {
    qclab::qgates::iSWAP< T >  iswap( 0 , 2 ) ;

    const T i(0,1) ;
    auto check = qclab::dense::zeros< T >( 8 ) ;
    check(0,0) = 1 ;
    check(4,1) = i ;
    check(2,2) = 1 ;
    check(6,3) = i ;
    check(1,4) = i ;
    check(5,5) = 1 ;
    check(3,6) = i ;
    check(7,7) = 1 ;

    // apply (3 qubits)
    auto mat3 = I3 ;
    iswap.apply( qclab::Side::Left , qclab::Op::NoTrans , 3 , mat3 ) ;
    test_qclab_qgates_iSWAP_check( mat3 , check ) ;
  }

}


/*
 * complex float
 */
TEST( qclab_qgates_iSWAP , complex_float ) {
  test_qclab_qgates_iSWAP< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_iSWAP , complex_double ) {
  test_qclab_qgates_iSWAP< std::complex< double > >() ;
}

