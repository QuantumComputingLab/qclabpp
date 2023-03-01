#include <gtest/gtest.h>
#include "qclab/dense/kron.hpp"
#include "qclab/qgates/PauliX.hpp"
#include "qclab/qgates/RotationY.hpp"

template <typename T>
void test_qclab_qgates_QGate1() {

  using R = qclab::real_t< T > ;
  using V = std::vector< T > ;
  using M = qclab::dense::SquareMatrix< T > ;

  const auto I1 = qclab::dense::eye< T >( 2 ) ;
  const auto I2 = qclab::dense::eye< T >( 4 ) ;
  const auto I3 = qclab::dense::eye< T >( 8 ) ;

  const V v1 = { 3 , 5 } ;
  const V v2 = { 3 , 5 , 2 , 7 } ;
  const V v3 = { 3 , 5 , 2 , 7 , 4 , 1 , 8 , 3 } ;

  // apply (NoTrans)
  {
    qclab::qgates::PauliX< T >  X( 0 ) ;

    // nbQubits = 1
    auto vec1 = v1 ;
    X.apply( qclab::Op::NoTrans , 1 , vec1 ) ;
    V check1 = { 5 , 3 } ;
    EXPECT_TRUE( vec1 == check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ; T* vec1_ = vec1.data() ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { X.apply_device( qclab::Op::NoTrans , 1 , vec1_ ) ; }
    EXPECT_TRUE( vec1 == check1 ) ;
  #endif

    // nbQubits = 2
    auto vec2 = v2 ;
    X.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 2 , 7 , 3 , 5 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { X.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    X.setQubit( 1 ) ;
    vec2 = v2 ;
    X.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    check2 = { 5 , 3 , 7 , 2 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { X.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    X.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 2 , 7 , 3 , 5 , 8 , 3 , 4 , 1 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { X.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    X.setQubit( 0 ) ;
    vec3 = v3 ;
    X.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 4 , 1 , 8 , 3 , 3 , 5 , 2 , 7 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { X.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    X.setQubit( 2 ) ;
    vec3 = v3 ;
    X.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 5 , 3 , 7 , 2 , 1 , 4 , 3 , 8 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { X.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif
  }

  // apply (Trans)
  {
    qclab::qgates::RotationY< T >  Y( 0 , R(0) , R(1) ) ;

    // nbQubits = 1
    auto vec1 = v1 ;
    Y.apply( qclab::Op::Trans , 1 , vec1 ) ;
    V check1 = { 5 , -3 } ;
    EXPECT_TRUE( vec1 == check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ; T* vec1_ = vec1.data() ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { Y.apply_device( qclab::Op::Trans , 1 , vec1_ ) ; }
    EXPECT_TRUE( vec1 == check1 ) ;
  #endif

    // nbQubits = 2
    auto vec2 = v2 ;
    Y.apply( qclab::Op::Trans , 2 , vec2 ) ;
    V check2 = { 2 , 7 , -3 , -5 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { Y.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    Y.setQubit( 1 ) ;
    vec2 = v2 ;
    Y.apply( qclab::Op::Trans , 2 , vec2 ) ;
    check2 = { 5 , -3 , 7 , -2 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { Y.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    Y.apply( qclab::Op::Trans , 3 , vec3 ) ;
    V check3 = { 2 , 7 , -3 , -5 , 8 , 3 , -4 , -1 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { Y.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    Y.setQubit( 0 ) ;
    vec3 = v3 ;
    Y.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { 4 , 1 , 8 , 3 , -3 , -5 , -2 , -7 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { Y.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    Y.setQubit( 2 ) ;
    vec3 = v3 ;
    Y.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { 5 , -3 , 7 , -2 , 1 , -4 , 3 , -8 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { Y.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif
  }


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

