#include <gtest/gtest.h>
#include "qclab/qgates/CNOT.hpp"

template <typename T, qclab::Side side>
void test_qclab_qgates_QControlledGate2() {

  using R = qclab::real_t< T > ;
  using V = std::vector< T > ;
  using M = qclab::dense::SquareMatrix< T > ;

  auto I1 = qclab::dense::eye< T >(  2 ) ;
  auto I2 = qclab::dense::eye< T >(  4 ) ;
  auto I3 = qclab::dense::eye< T >(  8 ) ;
  auto I4 = qclab::dense::eye< T >( 16 ) ;
  auto I5 = qclab::dense::eye< T >( 32 ) ;

  const V v2 = { 3 , 5 , 2 , 7 } ;
  const V v3 = { 3 , 5 , 2 , 7 , 4 , 1 , 8 , 3 } ;
  const V v4 = { 3 , 5 , 2 , 7 , 4 , 1 , 8 , 3 , 7 , 2 , 5 , 6 , 8 , 9 , 5 , 1};

  //
  // CNOTs |1> controlled
  //

  // control < target
  {
    qclab::qgates::CNOT< T >  cnot01( 0 , 1 ) ;
    qclab::qgates::CNOT< T >  cnot12( 1 , 2 ) ;
    qclab::qgates::CNOT< T >  cnot13( 1 , 3 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cnot01.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 3 , 5 , 7 , 2 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cnot01.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    auto mat2 = I2 ;
    cnot01.apply( side , qclab::Op::NoTrans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == cnot01.matrix() ) ;

    // nbQubits = 3
    auto vec3 = v3 ;
    cnot01.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 3 , 5 , 2 , 7 , 8 , 3 , 4 , 1 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cnot01.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif
    vec3 = v3 ;
    cnot12.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 7 , 2 , 4 , 1 , 3 , 8 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cnot12.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    auto mat3 = I3 ;
    cnot01.apply( side , qclab::Op::NoTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( cnot01.matrix() , I1 ) ) ;
    mat3 = I3 ;
    cnot12.apply( side , qclab::Op::NoTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( I1 , cnot12.matrix() ) ) ;

    // nbQubits = 4
    auto vec4 = v4 ;
    cnot12.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { 3 , 5 , 2 , 7 , 8 , 3 , 4 , 1 , 7 , 2 , 5 , 6 , 5 , 1 , 8 , 9};
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cnot12.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    auto mat4 = I4 ;
    cnot12.apply( side , qclab::Op::NoTrans , 4 , mat4 ) ;
    EXPECT_TRUE( mat4 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( cnot12.matrix() , I1 ) ) ) ;

    auto check = qclab::dense::zeros< T >( 8 ) ;
    check(0,0) = 1 ;
    check(1,1) = 1 ;
    check(2,2) = 1 ;
    check(3,3) = 1 ;
    check(5,4) = 1 ;
    check(4,5) = 1 ;
    check(7,6) = 1 ;
    check(6,7) = 1 ;

    // nbQubits = 5
    auto mat5 = I5 ;
    cnot13.apply( side , qclab::Op::NoTrans , 5 , mat5 ) ;
    EXPECT_TRUE( mat5 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( check , I1 ) ) ) ;
  }

  // control > target
  {
    qclab::qgates::CNOT< T >  cnot10( 1 , 0 ) ;
    qclab::qgates::CNOT< T >  cnot21( 2 , 1 ) ;
    qclab::qgates::CNOT< T >  cnot31( 3 , 1 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cnot10.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 3 , 7 , 2 , 5 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cnot10.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    auto mat2 = I2 ;
    cnot10.apply( side , qclab::Op::NoTrans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == cnot10.matrix() ) ;

    // nbQubits = 3
    auto vec3 = v3 ;
    cnot10.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 3 , 5 , 8 , 3 , 4 , 1 , 2 , 7 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cnot10.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif
    vec3 = v3 ;
    cnot21.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , 7 , 2 , 5 , 4 , 3 , 8 , 1 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cnot21.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    auto mat3 = I3 ;
    cnot10.apply( side , qclab::Op::NoTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( cnot10.matrix() , I1 ) ) ;
    mat3 = I3 ;
    cnot21.apply( side , qclab::Op::NoTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( I1 , cnot21.matrix() ) ) ;

    // nbQubits = 4
    auto vec4 = v4 ;
    cnot21.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { 3 , 5 , 8 , 3 , 4 , 1 , 2 , 7 , 7 , 2 , 5 , 1 , 8 , 9 , 5 , 6};
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cnot21.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    auto mat4 = I4 ;
    cnot21.apply( side , qclab::Op::NoTrans , 4 , mat4 ) ;
    EXPECT_TRUE( mat4 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( cnot21.matrix() , I1 ) ) ) ;

    auto check = qclab::dense::zeros< T >( 8 ) ;
    check(0,0) = 1 ;
    check(5,1) = 1 ;
    check(2,2) = 1 ;
    check(7,3) = 1 ;
    check(4,4) = 1 ;
    check(1,5) = 1 ;
    check(6,6) = 1 ;
    check(3,7) = 1 ;

    // nbQubits = 5
    auto mat5 = I5 ;
    cnot31.apply( side , qclab::Op::NoTrans , 5 , mat5 ) ;
    EXPECT_TRUE( mat5 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( check , I1 ) ) ) ;
  }


  //
  // CNOTs |0> controlled
  //

  // control < target
  {
    qclab::qgates::CNOT< T >  cnot01( 0 , 1 , 0 ) ;
    qclab::qgates::CNOT< T >  cnot12( 1 , 2 , 0 ) ;
    qclab::qgates::CNOT< T >  cnot13( 1 , 3 , 0 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cnot01.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 5 , 3 , 2 , 7 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cnot01.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    auto mat2 = I2 ;
    cnot01.apply( side , qclab::Op::NoTrans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == cnot01.matrix() ) ;

    // nbQubits = 3
    auto vec3 = v3 ;
    cnot01.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 2 , 7 , 3 , 5 , 4 , 1 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cnot01.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif
    vec3 = v3 ;
    cnot12.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 5 , 3 , 2 , 7 , 1 , 4 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cnot12.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    auto mat3 = I3 ;
    cnot01.apply( side , qclab::Op::NoTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( cnot01.matrix() , I1 ) ) ;
    mat3 = I3 ;
    cnot12.apply( side , qclab::Op::NoTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( I1 , cnot12.matrix() ) ) ;

    // nbQubits = 4
    auto vec4 = v4 ;
    cnot12.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { 2 , 7 , 3 , 5 , 4 , 1 , 8 , 3 , 5 , 6 , 7 , 2 , 8 , 9 , 5 , 1};
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cnot12.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    auto mat4 = I4 ;
    cnot12.apply( side , qclab::Op::NoTrans , 4 , mat4 ) ;
    EXPECT_TRUE( mat4 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( cnot12.matrix() , I1 ) ) ) ;

    auto check = qclab::dense::zeros< T >( 8 ) ;
    check(1,0) = 1 ;
    check(0,1) = 1 ;
    check(3,2) = 1 ;
    check(2,3) = 1 ;
    check(4,4) = 1 ;
    check(5,5) = 1 ;
    check(6,6) = 1 ;
    check(7,7) = 1 ;

    // nbQubits = 5
    auto mat5 = I5 ;
    cnot13.apply( side , qclab::Op::NoTrans , 5 , mat5 ) ;
    EXPECT_TRUE( mat5 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( check , I1 ) ) ) ;
  }

  // control > target
  {
    qclab::qgates::CNOT< T >  cnot10( 1 , 0 , 0 ) ;
    qclab::qgates::CNOT< T >  cnot21( 2 , 1 , 0 ) ;
    qclab::qgates::CNOT< T >  cnot31( 3 , 1 , 0 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cnot10.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 2 , 5 , 3 , 7 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cnot10.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    auto mat2 = I2 ;
    cnot10.apply( side , qclab::Op::NoTrans , 2 , mat2 ) ;
    EXPECT_TRUE( mat2 == cnot10.matrix() ) ;

    // nbQubits = 3
    auto vec3 = v3 ;
    cnot10.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 4 , 1 , 2 , 7 , 3 , 5 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cnot10.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif
    vec3 = v3 ;
    cnot21.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 2 , 5 , 3 , 7 , 8 , 1 , 4 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cnot21.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    auto mat3 = I3 ;
    cnot10.apply( side , qclab::Op::NoTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( cnot10.matrix() , I1 ) ) ;
    mat3 = I3 ;
    cnot21.apply( side , qclab::Op::NoTrans , 3 , mat3 ) ;
    EXPECT_TRUE( mat3 == qclab::dense::kron( I1 , cnot21.matrix() ) ) ;

    // nbQubits = 4
    auto vec4 = v4 ;
    cnot21.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { 4 , 1 , 2 , 7 , 3 , 5 , 8 , 3 , 8 , 9 , 5 , 6 , 7 , 2 , 5 , 1};
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cnot21.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    auto mat4 = I4 ;
    cnot21.apply( side , qclab::Op::NoTrans , 4 , mat4 ) ;
    EXPECT_TRUE( mat4 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( cnot21.matrix() , I1 ) ) ) ;

    auto check = qclab::dense::zeros< T >( 8 ) ;
    check(4,0) = 1 ;
    check(1,1) = 1 ;
    check(6,2) = 1 ;
    check(3,3) = 1 ;
    check(0,4) = 1 ;
    check(5,5) = 1 ;
    check(2,6) = 1 ;
    check(7,7) = 1 ;

    // nbQubits = 5
    auto mat5 = I5 ;
    cnot31.apply( side , qclab::Op::NoTrans , 5 , mat5 ) ;
    EXPECT_TRUE( mat5 == qclab::dense::kron( I1 ,
                           qclab::dense::kron( check , I1 ) ) ) ;
  }

}


/*
 * complex float
 */
TEST( qclab_qgates_QControlledGate2 , complex_float ) {
  test_qclab_qgates_QControlledGate2< std::complex< float > ,
                                      qclab::Side::Left >() ;
  test_qclab_qgates_QControlledGate2< std::complex< float > ,
                                      qclab::Side::Right >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_QControlledGate2 , complex_double ) {
  test_qclab_qgates_QControlledGate2< std::complex< double > ,
                                      qclab::Side::Left >() ;
  test_qclab_qgates_QControlledGate2< std::complex< double > ,
                                      qclab::Side::Right >() ;
}

