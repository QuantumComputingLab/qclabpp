#include <gtest/gtest.h>
#include "qclab/qgates/Phase90.hpp"

template <typename T>
void test_qclab_qgates_Phase90() {

  qclab::qgates::Phase90< T >  S ;

  EXPECT_EQ( S.nbQubits() , 1 ) ;   // nbQubits
  EXPECT_TRUE( S.fixed() ) ;        // fixed
  EXPECT_FALSE( S.controlled() ) ;  // controlled

  // qubit
  EXPECT_EQ( S.qubit() , 0 ) ;
  S.setQubit( 2 ) ;
  EXPECT_EQ( S.qubit() , 2 ) ;

  // qubits
  auto qubits = S.qubits() ;
  EXPECT_EQ( qubits.size() , 1 ) ;
  EXPECT_EQ( qubits[0] , 2 ) ;
  int qnew = 3 ;
  S.setQubits( &qnew ) ;
  EXPECT_EQ( S.qubit() , 3 ) ;

  // matrix
  EXPECT_EQ( std::real( S.matrix()(0,0) ) , 1.0 ) ;
  EXPECT_EQ( std::real( S.matrix()(1,0) ) , 0.0 ) ;
  EXPECT_EQ( std::real( S.matrix()(0,1) ) , 0.0 ) ;
  EXPECT_EQ( std::real( S.matrix()(1,1) ) , 0.0 ) ;
  EXPECT_EQ( std::imag( S.matrix()(0,0) ) , 0.0 ) ;
  EXPECT_EQ( std::imag( S.matrix()(1,0) ) , 0.0 ) ;
  EXPECT_EQ( std::imag( S.matrix()(0,1) ) , 0.0 ) ;
  EXPECT_EQ( std::imag( S.matrix()(1,1) ) , 1.0 ) ;

  // print
  S.print() ;

  // toQASM
  std::stringstream qasm ;
  EXPECT_EQ( S.toQASM( qasm ) , 0 ) ;
  EXPECT_EQ( qasm.str() , "s q[3];\n" ) ;
  std::cout << qasm.str() ;

  // operators == and !=
  qclab::qgates::Phase90< T > S2 ;
  EXPECT_TRUE( S == S2 ) ;
  EXPECT_FALSE( S != S2 ) ;

  // apply
  {
    using V = std::vector< T > ;
    const V v1 = { 1 , 2 } ;
    const V v2 = { 1 , 2 , 3 , 4 } ;
    const V v3 = { 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 } ;
    const V v4 = {  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  8 ,
                    9 , 10 , 11 , 12 , 13 , 14 , 15 , 16 } ;

    const T c  = T( 0 ,  1 ) ;
    const T cc = T( 0 , -1 ) ;

    qclab::qgates::Phase90< T >  P0( 0 ) ;
    qclab::qgates::Phase90< T >  P1( 1 ) ;
    qclab::qgates::Phase90< T >  P2( 2 ) ;
    qclab::qgates::Phase90< T >  P3( 3 ) ;

    // nbQubits = 1
    auto vec1 = v1 ;
    P0.apply( qclab::Op::NoTrans , 1 , vec1 ) ;
    V check1 = { 1 , T(2)*c } ;
    EXPECT_TRUE( vec1 == check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ; T* vec1_ = vec1.data() ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { P0.apply_device( qclab::Op::NoTrans , 1 , vec1_ ) ; }
    EXPECT_TRUE( vec1 == check1 ) ;
  #endif

    vec1 = v1 ;
    P0.apply( qclab::Op::Trans , 1 , vec1 ) ;
    EXPECT_TRUE( vec1 == check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { P0.apply_device( qclab::Op::Trans , 1 , vec1_ ) ; }
    EXPECT_TRUE( vec1 == check1 ) ;
  #endif

    vec1 = v1 ;
    P0.apply( qclab::Op::ConjTrans , 1 , vec1 ) ;
    check1 = { 1 , T(2)*cc } ;
    EXPECT_TRUE( vec1 == check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { P0.apply_device( qclab::Op::ConjTrans , 1 , vec1_ ) ; }
    EXPECT_TRUE( vec1 == check1 ) ;
  #endif

    // nbQubits = 2
    auto vec2 = v2 ;
    P0.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 1 , 2 , T(3)*c , T(4)*c } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { P0.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    P0.apply( qclab::Op::Trans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { P0.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    P0.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check2 = { 1 , 2 , T(3)*cc , T(4)*cc } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { P0.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    P1.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    check2 = { 1 , T(2)*c , 3 , T(4)*c } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { P1.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    P1.apply( qclab::Op::Trans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { P1.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    P1.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check2 = { 1 , T(2)*cc , 3 , T(4)*cc } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { P1.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    P0.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 1 , 2 , 3 , 4 , T(5)*c , T(6)*c , T(7)*c , T(8)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { P0.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    P0.apply( qclab::Op::Trans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { P0.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    P0.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 1 , 2 , 3 , 4 , T(5)*cc , T(6)*cc , T(7)*cc , T(8)*cc } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { P0.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    P1.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 1 , 2 , T(3)*c , T(4)*c , 5 , 6 , T(7)*c , T(8)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { P1.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    P1.apply( qclab::Op::Trans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { P1.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    P1.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 1 , 2 , T(3)*cc , T(4)*cc , 5 , 6 , T(7)*cc , T(8)*cc } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { P1.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    P2.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 1 , T(2)*c , 3 , T(4)*c , 5 , T(6)*c , 7 , T(8)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { P2.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    P2.apply( qclab::Op::Trans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { P2.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    P2.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 1 , T(2)*cc , 3 , T(4)*cc , 5 , T(6)*cc , 7 , T(8)*cc } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { P2.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    P0.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = {  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  8 ,
                 T( 9)*c , T(10)*c , T(11)*c , T(12)*c ,
                 T(13)*c , T(14)*c , T(15)*c , T(16)*c } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { P0.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    P0.apply( qclab::Op::Trans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { P0.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    P0.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = {  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  8 ,
               T( 9)*cc , T(10)*cc , T(11)*cc , T(12)*cc ,
               T(13)*cc , T(14)*cc , T(15)*cc , T(16)*cc } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { P0.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    P1.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = {  1 ,  2 ,  3 ,  4 , T( 5)*c , T( 6)*c , T( 7)*c , T( 8)*c ,
                9 , 10 , 11 , 12 , T(13)*c , T(14)*c , T(15)*c , T(16)*c } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { P1.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    P1.apply( qclab::Op::Trans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { P1.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    P1.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = {  1 ,  2 ,  3 ,  4 , T( 5)*cc , T( 6)*cc , T( 7)*cc , T( 8)*cc ,
                9 , 10 , 11 , 12 , T(13)*cc , T(14)*cc , T(15)*cc , T(16)*cc } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { P1.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    P2.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = {  1 ,  2 , T( 3)*c , T( 4)*c ,  5 ,  6 , T( 7)*c , T( 8)*c ,
                9 , 10 , T(11)*c , T(12)*c , 13 , 14 , T(15)*c , T(16)*c } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { P2.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    P2.apply( qclab::Op::Trans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { P2.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    P2.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = {  1 ,  2 , T( 3)*cc , T( 4)*cc ,  5 ,  6 , T( 7)*cc , T( 8)*cc ,
                9 , 10 , T(11)*cc , T(12)*cc , 13 , 14 , T(15)*cc , T(16)*cc } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { P2.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    P3.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = {  1 , T( 2)*c ,  3 , T( 4)*c ,  5 , T( 6)*c ,  7 , T( 8)*c ,
                9 , T(10)*c , 11 , T(12)*c , 13 , T(14)*c , 15 , T(16)*c } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { P3.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    P3.apply( qclab::Op::Trans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { P3.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    P3.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = {  1 , T( 2)*cc ,  3 , T( 4)*cc ,  5 , T( 6)*cc ,  7 , T( 8)*cc ,
                9 , T(10)*cc , 11 , T(12)*cc , 13 , T(14)*cc , 15 , T(16)*cc } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { P3.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif
  }

}


/*
 * complex float
 */
TEST( qclab_qgates_Phase90 , complex_float ) {
  test_qclab_qgates_Phase90< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_Phase90 , complex_double ) {
  test_qclab_qgates_Phase90< std::complex< double > >() ;
}

