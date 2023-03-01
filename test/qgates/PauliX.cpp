#include <gtest/gtest.h>
#include "qclab/qgates/PauliX.hpp"

template <typename T>
void test_qclab_qgates_PauliX() {

  qclab::qgates::PauliX< T >  X ;

  EXPECT_EQ( X.nbQubits() , 1 ) ;   // nbQubits
  EXPECT_TRUE( X.fixed() ) ;        // fixed
  EXPECT_FALSE( X.controlled() ) ;  // controlled

  // qubit
  EXPECT_EQ( X.qubit() , 0 ) ;
  X.setQubit( 2 ) ;
  EXPECT_EQ( X.qubit() , 2 ) ;

  // qubits
  auto qubits = X.qubits() ;
  EXPECT_EQ( qubits.size() , 1 ) ;
  EXPECT_EQ( qubits[0] , 2 ) ;
  int qnew = 3 ;
  X.setQubits( &qnew ) ;
  EXPECT_EQ( X.qubit() , 3 ) ;

  // matrix
  EXPECT_EQ( X.matrix()(0,0) , T(0.0) ) ;
  EXPECT_EQ( X.matrix()(1,0) , T(1.0) ) ;
  EXPECT_EQ( X.matrix()(0,1) , T(1.0) ) ;
  EXPECT_EQ( X.matrix()(1,1) , T(0.0) ) ;

  // print
  X.print() ;

  // toQASM
  std::stringstream qasm ;
  EXPECT_EQ( X.toQASM( qasm ) , 0 ) ;
  EXPECT_EQ( qasm.str() , "x q[3];\n" ) ;
  std::cout << qasm.str() ;

  // operators == and !=
  qclab::qgates::PauliX< T > X2 ;
  EXPECT_TRUE( X == X2 ) ;
  EXPECT_FALSE( X != X2 ) ;

  // apply
  {
    using V = std::vector< T > ;
    const V v1 = { 0 , 1 } ;
    const V v2 = { 0 , 1 , 2 , 3 } ;
    const V v3 = { 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 } ;
    const V v4 = {  0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,
                    8 ,  9 , 10 , 11 , 12 , 13 , 14 , 15 } ;

    qclab::qgates::PauliX< T >  X0( 0 ) ;
    qclab::qgates::PauliX< T >  X1( 1 ) ;
    qclab::qgates::PauliX< T >  X2( 2 ) ;
    qclab::qgates::PauliX< T >  X3( 3 ) ;

    // nbQubits = 1
    auto vec1 = v1 ;
    X0.apply( qclab::Op::NoTrans , 1 , vec1 ) ;
    V check1 = { 1 , 0 } ;
    EXPECT_TRUE( vec1 == check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ; T* vec1_ = vec1.data() ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { X0.apply_device( qclab::Op::NoTrans , 1 , vec1_ ) ; }
    EXPECT_TRUE( vec1 == check1 ) ;
  #endif

    vec1 = v1 ;
    X0.apply( qclab::Op::Trans , 1 , vec1 ) ;
    EXPECT_TRUE( vec1 == check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { X0.apply_device( qclab::Op::Trans , 1 , vec1_ ) ; }
    EXPECT_TRUE( vec1 == check1 ) ;
  #endif

    vec1 = v1 ;
    X0.apply( qclab::Op::ConjTrans , 1 , vec1 ) ;
    EXPECT_TRUE( vec1 == check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { X0.apply_device( qclab::Op::ConjTrans , 1 , vec1_ ) ; }
    EXPECT_TRUE( vec1 == check1 ) ;
  #endif

    // nbQubits = 2
    auto vec2 = v2 ;
    X0.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 2 , 3 , 0 , 1 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { X0.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    X0.apply( qclab::Op::Trans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { X0.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    X0.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { X0.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    X1.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    check2 = { 1 , 0 , 3 , 2 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { X1.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    X1.apply( qclab::Op::Trans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { X1.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    X1.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { X1.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    X0.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 4 , 5 , 6 , 7 , 0 , 1 , 2 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { X0.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    X0.apply( qclab::Op::Trans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { X0.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    X0.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { X0.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    X1.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 2 , 3 , 0 , 1 , 6 , 7 , 4 , 5 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { X1.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    X1.apply( qclab::Op::Trans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { X1.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    X1.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { X1.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    X2.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 1 , 0 , 3 , 2 , 5 , 4 , 7 , 6 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { X2.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    X2.apply( qclab::Op::Trans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { X2.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    X2.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { X2.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    X0.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = {  8 ,  9 , 10 , 11 , 12 , 13 , 14 , 15 ,
                  0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { X0.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    X0.apply( qclab::Op::Trans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { X0.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    X0.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { X0.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    X1.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = {  4 ,  5 ,  6 ,  7 ,  0 ,  1 ,  2 ,  3 ,
               12 , 13 , 14 , 15 ,  8 ,  9 , 10 , 11 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { X1.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    X1.apply( qclab::Op::Trans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { X1.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    X1.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { X1.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    X2.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = {  2 ,  3 ,  0 ,  1 ,  6 ,  7 ,  4 ,  5 ,
               10 , 11 ,  8 ,  9 , 14 , 15 , 12 , 13 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { X2.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    X2.apply( qclab::Op::Trans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { X2.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    X2.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { X2.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    X3.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = {  1 ,  0 ,  3 ,  2 ,  5 ,  4 ,  7 ,  6 ,
                9 ,  8 , 11 , 10 , 13 , 12 , 15 , 14 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { X3.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    X3.apply( qclab::Op::Trans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { X3.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    X3.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { X3.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif
  }

}


/*
 * float
 */
TEST( qclab_qgates_PauliX , float ) {
  test_qclab_qgates_PauliX< float >() ;
}

/*
 * double
 */
TEST( qclab_qgates_PauliX , double ) {
  test_qclab_qgates_PauliX< double >() ;
}


/*
 * complex float
 */
TEST( qclab_qgates_PauliX , complex_float ) {
  test_qclab_qgates_PauliX< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_PauliX , complex_double ) {
  test_qclab_qgates_PauliX< std::complex< double > >() ;
}

