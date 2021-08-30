#include <gtest/gtest.h>
#include "qclab/qgates/RotationXX.hpp"
#include "qclab/qgates/RotationYY.hpp"
#include "qclab/dense/kron.hpp"

template <typename T>
void test_qclab_qgates_QGate2() {

  using R = qclab::real_t< T > ;
  using V = std::vector< T > ;
  using M = qclab::dense::SquareMatrix< T > ;

  const auto I1 = qclab::dense::eye< T >(  2 ) ;
  const auto I2 = qclab::dense::eye< T >(  4 ) ;
  const auto I3 = qclab::dense::eye< T >(  8 ) ;
  const auto I4 = qclab::dense::eye< T >( 16 ) ;

  const V v2 = { 3 , 5 , 2 , 7 } ;
  const V v3 = { 3 , 5 , 2 , 7 , 4 , 1 , 8 , 3 } ;
  const V v4 = { 3 , 5 , 2 , 7 , 4 , 1 , 8 , 3 , 7 , 2 , 5 , 6 , 8 , 9 , 5 , 1};

  // apply (NoTrans)
  {
    qclab::qgates::RotationXX< T >  XX( 0 , 1 , R(0) , R(1) ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    XX.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { T(0,-7) , T(0,-2) , T(0,-5) , T(0,-3) } ;
    EXPECT_TRUE( vec2 == check2 ) ;

    // nbQubits = 3
    auto vec3 = v3 ;
    XX.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { T(0,-8) , T(0,-3) , T(0,-4) , T(0,-1) ,
                 T(0,-2) , T(0,-7) , T(0,-3) , T(0,-5) } ;
    EXPECT_TRUE( vec3 == check3 ) ;

    int qubits[2] = { 1 , 2 } ;
    XX.setQubits( &qubits[0] ) ;
    vec3 = v3 ;
    XX.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T(0,-7) , T(0,-2) , T(0,-5) , T(0,-3) ,
               T(0,-3) , T(0,-8) , T(0,-1) , T(0,-4) } ;
    EXPECT_TRUE( vec3 == check3 ) ;

    // nbQubits = 4
    auto vec4 = v4 ;
    XX.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { T(0,-8) , T(0,-3) , T(0,-4) , T(0,-1) ,
                 T(0,-2) , T(0,-7) , T(0,-3) , T(0,-5) ,
                 T(0,-5) , T(0,-1) , T(0,-8) , T(0,-9) ,
                 T(0,-5) , T(0,-6) , T(0,-7) , T(0,-2) } ;
    EXPECT_TRUE( vec4 == check4 ) ;

    qubits[0] = 0 ;
    qubits[1] = 1 ;
    XX.setQubits( &qubits[0] ) ;
    vec4 = v4 ;
    XX.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = { T(0,-8) , T(0,-9) , T(0,-5) , T(0,-1) ,
               T(0,-7) , T(0,-2) , T(0,-5) , T(0,-6) ,
               T(0,-4) , T(0,-1) , T(0,-8) , T(0,-3) ,
               T(0,-3) , T(0,-5) , T(0,-2) , T(0,-7) } ;
    EXPECT_TRUE( vec4 == check4 ) ;

    qubits[0] = 2 ;
    qubits[1] = 3 ;
    XX.setQubits( &qubits[0] ) ;
    vec4 = v4 ;
    XX.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = { T(0,-7) , T(0,-2) , T(0,-5) , T(0,-3) ,
               T(0,-3) , T(0,-8) , T(0,-1) , T(0,-4) ,
               T(0,-6) , T(0,-5) , T(0,-2) , T(0,-7) ,
               T(0,-1) , T(0,-5) , T(0,-9) , T(0,-8) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  }

  // apply (Trans)
  {
    qclab::qgates::RotationYY< T >  YY( 0 , 1 , R(0) , R(1) ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    YY.apply( qclab::Op::Trans , 2 , vec2 ) ;
    V check2 = { T(0,7) , T(0,-2) , T(0,-5) , T(0,3) } ;
    EXPECT_TRUE( vec2 == check2 ) ;

    // nbQubits = 3
    auto vec3 = v3 ;
    YY.apply( qclab::Op::Trans , 3 , vec3 ) ;
    V check3 = { T(0, 8) , T(0, 3) , T(0,-4) , T(0,-1) ,
                 T(0,-2) , T(0,-7) , T(0, 3) , T(0, 5) } ;
    EXPECT_TRUE( vec3 == check3 ) ;

    int qubits[2] = { 1 , 2 } ;
    YY.setQubits( &qubits[0] ) ;
    vec3 = v3 ;
    YY.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { T(0, 7) , T(0,-2) , T(0,-5) , T(0, 3) ,
               T(0, 3) , T(0,-8) , T(0,-1) , T(0, 4) } ;
    EXPECT_TRUE( vec3 == check3 ) ;

    // nbQubits = 4
    auto vec4 = v4 ;
    YY.apply( qclab::Op::Trans , 4 , vec4 ) ;
    V check4 = { T(0, 8) , T(0, 3) , T(0,-4) , T(0,-1) ,
                 T(0,-2) , T(0,-7) , T(0, 3) , T(0, 5) ,
                 T(0, 5) , T(0, 1) , T(0,-8) , T(0,-9) ,
                 T(0,-5) , T(0,-6) , T(0, 7) , T(0, 2) } ;
    EXPECT_TRUE( vec4 == check4 ) ;

    qubits[0] = 0 ;
    qubits[1] = 1 ;
    YY.setQubits( &qubits[0] ) ;
    vec4 = v4 ;
    YY.apply( qclab::Op::Trans , 4 , vec4 ) ;
    check4 = { T(0, 8) , T(0, 9) , T(0, 5) , T(0, 1) ,
               T(0,-7) , T(0,-2) , T(0,-5) , T(0,-6) ,
               T(0,-4) , T(0,-1) , T(0,-8) , T(0,-3) ,
               T(0, 3) , T(0, 5) , T(0, 2) , T(0, 7) } ;
    EXPECT_TRUE( vec4 == check4 ) ;

    qubits[0] = 2 ;
    qubits[1] = 3 ;
    YY.setQubits( &qubits[0] ) ;
    vec4 = v4 ;
    YY.apply( qclab::Op::Trans , 4 , vec4 ) ;
    check4 = { T(0, 7) , T(0,-2) , T(0,-5) , T(0, 3) ,
               T(0, 3) , T(0,-8) , T(0,-1) , T(0, 4) ,
               T(0, 6) , T(0,-5) , T(0,-2) , T(0, 7) ,
               T(0, 1) , T(0,-5) , T(0,-9) , T(0, 8) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  }

  // apply (ConjTrans)
  {
    qclab::qgates::RotationYY< T >  YY( 0 , 1 , R(0) , R(1) ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    YY.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    V check2 = { T(0,-7) , T(0,2) , T(0,5) , T(0,-3) } ;
    EXPECT_TRUE( vec2 == check2 ) ;

    // nbQubits = 3
    auto vec3 = v3 ;
    YY.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    V check3 = { T(0,-8) , T(0,-3) , T(0, 4) , T(0, 1) ,
                 T(0, 2) , T(0, 7) , T(0,-3) , T(0,-5) } ;
    EXPECT_TRUE( vec3 == check3 ) ;

    int qubits[2] = { 1 , 2 } ;
    YY.setQubits( &qubits[0] ) ;
    vec3 = v3 ;
    YY.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { T(0,-7) , T(0, 2) , T(0, 5) , T(0,-3) ,
               T(0,-3) , T(0, 8) , T(0, 1) , T(0,-4) } ;
    EXPECT_TRUE( vec3 == check3 ) ;

    // nbQubits = 4
    auto vec4 = v4 ;
    YY.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    V check4 = { T(0,-8) , T(0,-3) , T(0, 4) , T(0, 1) ,
                 T(0, 2) , T(0, 7) , T(0,-3) , T(0,-5) ,
                 T(0,-5) , T(0,-1) , T(0, 8) , T(0, 9) ,
                 T(0, 5) , T(0, 6) , T(0,-7) , T(0,-2) } ;
    EXPECT_TRUE( vec4 == check4 ) ;

    qubits[0] = 0 ;
    qubits[1] = 1 ;
    YY.setQubits( &qubits[0] ) ;
    vec4 = v4 ;
    YY.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { T(0,-8) , T(0,-9) , T(0,-5) , T(0,-1) ,
               T(0, 7) , T(0, 2) , T(0, 5) , T(0, 6) ,
               T(0, 4) , T(0, 1) , T(0, 8) , T(0, 3) ,
               T(0,-3) , T(0,-5) , T(0,-2) , T(0,-7) } ;
    EXPECT_TRUE( vec4 == check4 ) ;

    qubits[0] = 2 ;
    qubits[1] = 3 ;
    YY.setQubits( &qubits[0] ) ;
    vec4 = v4 ;
    YY.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { T(0,-7) , T(0, 2) , T(0, 5) , T(0,-3) ,
               T(0,-3) , T(0, 8) , T(0, 1) , T(0,-4) ,
               T(0,-6) , T(0, 5) , T(0, 2) , T(0,-7) ,
               T(0,-1) , T(0, 5) , T(0, 9) , T(0,-8) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  }


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

