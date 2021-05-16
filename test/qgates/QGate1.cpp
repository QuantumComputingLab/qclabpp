#include <gtest/gtest.h>
#include "qgates/PauliX.hpp"
#include "qgates/RotationY.hpp"

template <typename T>
void test_qclab_qgates_QGate1() {

  using R = qclab::real_t< T > ;
  using M = qclab::dense::SquareMatrix< T > ;

  auto I1 = qclab::dense::eye< T >( 2 ) ;
  auto I2 = qclab::dense::eye< T >( 4 ) ;
  auto I3 = qclab::dense::eye< T >( 8 ) ;

  // apply (Left + NoTrans)
  {
    qclab::qgates::PauliX< T >  X( 0 ) ;

    // nbQubits = 1
    auto mat1 = I1 ;
    X.apply( qclab::Side::Left , qclab::Op::NoTrans , 1 , mat1 ) ;
    EXPECT_TRUE( mat1 == X.matrix() ) ;
    mat1 = M( 2 , 3 ,
              5 , 7 ) ;
    X.apply( qclab::Side::Left , qclab::Op::NoTrans , 1 , mat1 ) ;
    M check1( 3 , 2 ,
              7 , 5 ) ;
    EXPECT_TRUE( mat1 == check1 ) ;

    // nbQubits = 2
    auto mat2 = I2 ;
    X.apply( qclab::Side::Left , qclab::Op::NoTrans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == qclab::dense::kron( X.matrix() , I1 ) ) ;

    X.setQubit( 1 ) ;
    mat2 = I2 ;
    X.apply( qclab::Side::Left , qclab::Op::NoTrans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == qclab::dense::kron( I1 , X.matrix() ) ) ;

    // nbQubits = 3
    auto mat3 = I3 ;
    X.apply( qclab::Side::Left , qclab::Op::NoTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( X.matrix() , I1 ) ) ) ;
  }

  // apply (Right + NoTrans)
  {
    qclab::qgates::PauliX< T >  X( 0 ) ;

    // nbQubits = 1
    auto mat1 = I1 ;
    X.apply( qclab::Side::Right , qclab::Op::NoTrans , 1 , mat1 ) ;
    EXPECT_TRUE( mat1 == X.matrix() ) ;
    mat1 = M( 2 , 3 ,
              5 , 7 ) ;
    X.apply( qclab::Side::Right , qclab::Op::NoTrans , 1 , mat1 ) ;
    M check1( 5 , 7 ,
              2 , 3 ) ;
    EXPECT_TRUE( mat1 == check1 ) ;

    // nbQubits = 2
    auto mat2 = I2 ;
    X.apply( qclab::Side::Right , qclab::Op::NoTrans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == qclab::dense::kron( X.matrix() , I1 ) ) ;

    X.setQubit( 1 ) ;
    mat2 = I2 ;
    X.apply( qclab::Side::Right , qclab::Op::NoTrans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == qclab::dense::kron( I1 , X.matrix() ) ) ;

    // nbQubits = 3
    auto mat3 = I3 ;
    X.apply( qclab::Side::Right , qclab::Op::NoTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( X.matrix() , I1 ) ) ) ;
  }

  // apply (Left + Trans)
  {
    qclab::qgates::RotationY< T >  Y( 0 , R(0) , R(1) ) ;
    M Ytrans = qclab::dense::trans( Y.matrix() ) ;

    // nbQubits = 1
    auto mat1 = I1 ;
    Y.apply( qclab::Side::Left , qclab::Op::Trans , 1 , mat1 ) ;
    EXPECT_TRUE( mat1 == Ytrans ) ;
    mat1 = M( 2 , 3 ,
              5 , 7 ) ;
    Y.apply( qclab::Side::Left , qclab::Op::Trans , 1 , mat1 ) ;
    M check1( -3 , 2 ,
              -7 , 5 ) ;
    EXPECT_TRUE( mat1 == check1 ) ;

    // nbQubits = 2
    auto mat2 = I2 ;
    Y.apply( qclab::Side::Left , qclab::Op::Trans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == qclab::dense::kron( Ytrans , I1 ) ) ;

    Y.setQubit( 1 ) ;
    mat2 = I2 ;
    Y.apply( qclab::Side::Left , qclab::Op::Trans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == qclab::dense::kron( I1 , Ytrans ) ) ;

    // nbQubits = 3
    auto mat3 = I3 ;
    Y.apply( qclab::Side::Left , qclab::Op::Trans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( Ytrans , I1 ) ) ) ;
  }

  // apply (Right + Trans)
  {
    qclab::qgates::RotationY< T >  Y( 0 , R(0) , R(1) ) ;
    M Ytrans = qclab::dense::trans( Y.matrix() ) ;

    // nbQubits = 1
    auto mat1 = I1 ;
    Y.apply( qclab::Side::Right , qclab::Op::Trans , 1 , mat1 ) ;
    EXPECT_TRUE( mat1 == Ytrans ) ;
    mat1 = M( 2 , 3 ,
              5 , 7 ) ;
    Y.apply( qclab::Side::Right , qclab::Op::Trans , 1 , mat1 ) ;
    M check1(  5 ,  7 ,
              -2 , -3 ) ;
    EXPECT_TRUE( mat1 == check1 ) ;

    // nbQubits = 2
    auto mat2 = I2 ;
    Y.apply( qclab::Side::Right , qclab::Op::Trans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == qclab::dense::kron( Ytrans , I1 ) ) ;

    Y.setQubit( 1 ) ;
    mat2 = I2 ;
    Y.apply( qclab::Side::Right , qclab::Op::Trans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == qclab::dense::kron( I1 , Ytrans ) ) ;

    // nbQubits = 3
    auto mat3 = I3 ;
    Y.apply( qclab::Side::Right , qclab::Op::Trans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( Ytrans , I1 ) ) ) ;
  }

}


/*
 * float
 */
TEST( qclab_qgates_QGate1 , float ) {
  test_qclab_qgates_QGate1< float >() ;
}

/*
 * double
 */
TEST( qclab_qgates_QGate1 , double ) {
  test_qclab_qgates_QGate1< double >() ;
}


/*
 * complex float
 */
TEST( qclab_qgates_QGate1 , complex_float ) {
  test_qclab_qgates_QGate1< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_QGate1 , complex_double ) {
  test_qclab_qgates_QGate1< std::complex< double > >() ;
}

