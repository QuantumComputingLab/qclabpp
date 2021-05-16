#include <gtest/gtest.h>
#include "qgates/RotationXX.hpp"
#include "qgates/RotationYY.hpp"

template <typename T>
void test_qclab_qgates_QGate2() {

  using R = qclab::real_t< T > ;
  using M = qclab::dense::SquareMatrix< T > ;

  auto I1 = qclab::dense::eye< T >(  2 ) ;
  auto I2 = qclab::dense::eye< T >(  4 ) ;
  auto I3 = qclab::dense::eye< T >(  8 ) ;
  auto I4 = qclab::dense::eye< T >( 16 ) ;

  // apply (Left + NoTrans)
  {
    qclab::qgates::RotationXX< T >  XX( 0 , 1 , R(0) , R(1) ) ;

    // nbQubits = 2
    auto mat2 = I2 ;
    XX.apply( qclab::Side::Left , qclab::Op::NoTrans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == XX.matrix() ) ;

    // nbQubits = 3
    auto mat3 = I3 ;
    XX.apply( qclab::Side::Left , qclab::Op::NoTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( XX.matrix() , I1 ) ) ;

    int qubits[2] = { 1 , 2 } ;
    XX.setQubits( &qubits[0] ) ;
    mat3 = I3 ;
    XX.apply( qclab::Side::Left , qclab::Op::NoTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( I1 , XX.matrix() ) ) ;

    // nbQubits = 4
    auto mat4 = I4 ;
    XX.apply( qclab::Side::Left , qclab::Op::NoTrans , 4 , mat4 ) ;
    EXPECT_TRUE( mat4 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( XX.matrix() , I1 ) ) ) ;
  }

  // apply (Right + NoTrans)
  {
    qclab::qgates::RotationXX< T >  XX( 0 , 1 , R(0) , R(1) ) ;

    // nbQubits = 2
    auto mat2 = I2 ;
    XX.apply( qclab::Side::Right , qclab::Op::NoTrans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == XX.matrix() ) ;

    // nbQubits = 3
    auto mat3 = I3 ;
    XX.apply( qclab::Side::Right , qclab::Op::NoTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( XX.matrix() , I1 ) ) ;

    int qubits[2] = { 1 , 2 } ;
    XX.setQubits( &qubits[0] ) ;
    mat3 = I3 ;
    XX.apply( qclab::Side::Right , qclab::Op::NoTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( I1 , XX.matrix() ) ) ;

    // nbQubits = 4
    auto mat4 = I4 ;
    XX.apply( qclab::Side::Right , qclab::Op::NoTrans , 4 , mat4 ) ;
    EXPECT_TRUE( mat4 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( XX.matrix() , I1 ) ) ) ;
  }

  // apply (Left + Trans)
  {
    qclab::qgates::RotationYY< T >  YY( 0 , 1 , R(0) , R(1) ) ;
    M YYtrans = qclab::dense::trans( YY.matrix() ) ;

    // nbQubits = 2
    auto mat2 = I2 ;
    YY.apply( qclab::Side::Left , qclab::Op::Trans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == YYtrans ) ;

    // nbQubits = 3
    auto mat3 = I3 ;
    YY.apply( qclab::Side::Left , qclab::Op::Trans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( YYtrans , I1 ) ) ;

    int qubits[2] = { 1 , 2 } ;
    YY.setQubits( &qubits[0] ) ;
    mat3 = I3 ;
    YY.apply( qclab::Side::Left , qclab::Op::Trans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( I1 , YYtrans ) ) ;

    // nbQubits = 4
    auto mat4 = I4 ;
    YY.apply( qclab::Side::Left , qclab::Op::Trans , 4 , mat4 ) ;
    EXPECT_TRUE( mat4 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( YYtrans , I1 ) ) ) ;
  }

  // apply (Right + Trans)
  {
    qclab::qgates::RotationYY< T >  YY( 0 , 1 , R(0) , R(1) ) ;
    M YYtrans = qclab::dense::trans( YY.matrix() ) ;

    // nbQubits = 2
    auto mat2 = I2 ;
    YY.apply( qclab::Side::Right , qclab::Op::Trans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == YYtrans ) ;

    // nbQubits = 3
    auto mat3 = I3 ;
    YY.apply( qclab::Side::Right , qclab::Op::Trans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( YYtrans , I1 ) ) ;

    int qubits[2] = { 1 , 2 } ;
    YY.setQubits( &qubits[0] ) ;
    mat3 = I3 ;
    YY.apply( qclab::Side::Right , qclab::Op::Trans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( I1 , YYtrans ) ) ;

    // nbQubits = 4
    auto mat4 = I4 ;
    YY.apply( qclab::Side::Right , qclab::Op::Trans , 4 , mat4 ) ;
    EXPECT_TRUE( mat4 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( YYtrans , I1 ) ) ) ;
  }

  // apply (Left + ConjTrans)
  {
    qclab::qgates::RotationYY< T >  YY( 0 , 1 , R(0) , R(1) ) ;
    M YYconjTrans = qclab::dense::conjTrans( YY.matrix() ) ;

    // nbQubits = 2
    auto mat2 = I2 ;
    YY.apply( qclab::Side::Left , qclab::Op::ConjTrans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == YYconjTrans ) ;

    // nbQubits = 3
    auto mat3 = I3 ;
    YY.apply( qclab::Side::Left , qclab::Op::ConjTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( YYconjTrans , I1 ) ) ;

    int qubits[2] = { 1 , 2 } ;
    YY.setQubits( &qubits[0] ) ;
    mat3 = I3 ;
    YY.apply( qclab::Side::Left , qclab::Op::ConjTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( I1 , YYconjTrans ) ) ;

    // nbQubits = 4
    auto mat4 = I4 ;
    YY.apply( qclab::Side::Left , qclab::Op::ConjTrans , 4 , mat4 ) ;
    EXPECT_TRUE( mat4 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( YYconjTrans , I1 ) ) ) ;
  }

  // apply (Right + ConjTrans)
  {
    qclab::qgates::RotationYY< T >  YY( 0 , 1 , R(0) , R(1) ) ;
    M YYconjConjTrans = qclab::dense::conjTrans( YY.matrix() ) ;

    // nbQubits = 2
    auto mat2 = I2 ;
    YY.apply( qclab::Side::Right , qclab::Op::ConjTrans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == YYconjConjTrans ) ;

    // nbQubits = 3
    auto mat3 = I3 ;
    YY.apply( qclab::Side::Right , qclab::Op::ConjTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( YYconjConjTrans , I1 ) ) ;

    int qubits[2] = { 1 , 2 } ;
    YY.setQubits( &qubits[0] ) ;
    mat3 = I3 ;
    YY.apply( qclab::Side::Right , qclab::Op::ConjTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( I1 , YYconjConjTrans ) ) ;

    // nbQubits = 4
    auto mat4 = I4 ;
    YY.apply( qclab::Side::Right , qclab::Op::ConjTrans , 4 , mat4 ) ;
    EXPECT_TRUE( mat4 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( YYconjConjTrans , I1 ) ) ) ;
  }

}


/*
 * complex float
 */
TEST( qclab_qgates_QGate2 , complex_float ) {
  test_qclab_qgates_QGate2< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_QGate2 , complex_double ) {
  test_qclab_qgates_QGate2< std::complex< double > >() ;
}

