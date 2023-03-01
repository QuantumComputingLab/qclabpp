#include <gtest/gtest.h>
#include "qclab/qgates/PauliY.hpp"

template <typename T>
void test_qclab_qgates_PauliY() {

  qclab::qgates::PauliY< T >  Y ;

  EXPECT_EQ( Y.nbQubits() , 1 ) ;   // nbQubits
  EXPECT_TRUE( Y.fixed() ) ;        // fixed
  EXPECT_FALSE( Y.controlled() ) ;  // controlled

  // qubit
  EXPECT_EQ( Y.qubit() , 0 ) ;
  Y.setQubit( 2 ) ;
  EXPECT_EQ( Y.qubit() , 2 ) ;

  // qubits
  auto qubits = Y.qubits() ;
  EXPECT_EQ( qubits.size() , 1 ) ;
  EXPECT_EQ( qubits[0] , 2 ) ;
  int qnew = 3 ;
  Y.setQubits( &qnew ) ;
  EXPECT_EQ( Y.qubit() , 3 ) ;

  // matrix
  EXPECT_EQ( std::real( Y.matrix()(0,0) ) ,  0.0 ) ;
  EXPECT_EQ( std::real( Y.matrix()(1,0) ) ,  0.0 ) ;
  EXPECT_EQ( std::real( Y.matrix()(0,1) ) ,  0.0 ) ;
  EXPECT_EQ( std::real( Y.matrix()(1,1) ) ,  0.0 ) ;
  EXPECT_EQ( std::imag( Y.matrix()(0,0) ) ,  0.0 ) ;
  EXPECT_EQ( std::imag( Y.matrix()(1,0) ) ,  1.0 ) ;
  EXPECT_EQ( std::imag( Y.matrix()(0,1) ) , -1.0 ) ;
  EXPECT_EQ( std::imag( Y.matrix()(1,1) ) ,  0.0 ) ;

  // print
  Y.print() ;

  // toQASM
  std::stringstream qasm ;
  EXPECT_EQ( Y.toQASM( qasm ) , 0 ) ;
  EXPECT_EQ( qasm.str() , "y q[3];\n" ) ;
  std::cout << qasm.str() ;

  // operators == and !=
  qclab::qgates::PauliY< T > Y2 ;
  EXPECT_TRUE( Y == Y2 ) ;
  EXPECT_FALSE( Y != Y2 ) ;

  // apply
  {
    using V = std::vector< T > ;
    const V v1 = { T(1,1) , T(2,2) } ;
    const V v2 = { T(1,1) , T(2,2) , T(3,3) , T(4,4) } ;
    const V v3 = { T(1,1) , T(2,2) , T(3,3) , T(4,4) ,
                   T(5,5) , T(6,6) , T(7,7) , T(8,8) } ;
    const V v4 = { T( 1, 1) , T( 2, 2) , T( 3, 3) , T( 4, 4) ,
                   T( 5, 5) , T( 6, 6) , T( 7, 7) , T( 8, 8) ,
                   T( 9, 9) , T(10,10) , T(11,11) , T(12,12) ,
                   T(13,13) , T(14,14) , T(15,15) , T(16,16) } ;

    qclab::qgates::PauliY< T >  Y0( 0 ) ;
    qclab::qgates::PauliY< T >  Y1( 1 ) ;
    qclab::qgates::PauliY< T >  Y2( 2 ) ;
    qclab::qgates::PauliY< T >  Y3( 3 ) ;

    // nbQubits = 1
    auto vec1 = v1 ;
    Y0.apply( qclab::Op::NoTrans , 1 , vec1 ) ;
    V check1 = { T( 2,-2) , T(-1, 1) } ;
    EXPECT_TRUE( vec1 == check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ; T* vec1_ = vec1.data() ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { Y0.apply_device( qclab::Op::NoTrans , 1 , vec1_ ) ; }
    EXPECT_TRUE( vec1 == check1 ) ;
  #endif

    vec1 = v1 ;
    Y0.apply( qclab::Op::ConjTrans , 1 , vec1 ) ;
    EXPECT_TRUE( vec1 == check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ;
    #pragma omp target data map(tofrom:vec1_[0:2])
   { Y0.apply_device( qclab::Op::ConjTrans , 1 , vec1_ ) ; }
    EXPECT_TRUE( vec1 == check1 ) ;
  #endif

    vec1 = v1 ;
    Y0.apply( qclab::Op::Trans , 1 , vec1 ) ;
    check1 = { T(-2, 2) , T( 1,-1) } ;
    EXPECT_TRUE( vec1 == check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { Y0.apply_device( qclab::Op::Trans , 1 , vec1_ ) ; }
    EXPECT_TRUE( vec1 == check1 ) ;
  #endif

    // nbQubits = 2
    auto vec2 = v2 ;
    Y0.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { T(3,-3) , T(4,-4) , T(-1,1) , T(-2,2) } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { Y0.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    Y0.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
   { Y0.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    Y0.apply( qclab::Op::Trans , 2 , vec2 ) ;
    check2 = { T(-3,3) , T(-4,4) , T(1,-1) , T(2,-2) } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { Y0.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    Y1.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    check2 = { T(2,-2) , T(-1,1) , T(4,-4) , T(-3,3) } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { Y1.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    Y1.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
   { Y1.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    Y1.apply( qclab::Op::Trans , 2 , vec2 ) ;
    check2 = { T(-2,2) , T(1,-1) , T(-4,4) , T(3,-3) } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { Y1.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    Y0.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { T( 5,-5) , T( 6,-6) , T( 7,-7) , T( 8,-8) ,
                 T(-1, 1) , T(-2, 2) , T(-3, 3) , T(-4, 4) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { Y0.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    Y0.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
   { Y0.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    Y0.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { T(-5, 5) , T(-6, 6) , T(-7, 7) , T(-8, 8) ,
               T( 1,-1) , T( 2,-2) , T( 3,-3) , T( 4,-4) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { Y0.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    Y1.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T( 3,-3) , T( 4,-4) , T(-1, 1) , T(-2, 2) ,
               T( 7,-7) , T( 8,-8) , T(-5, 5) , T(-6, 6) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { Y1.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    Y1.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
   { Y1.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    Y1.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { T(-3, 3) , T(-4, 4) , T( 1,-1) , T( 2,-2) ,
               T(-7, 7) , T(-8, 8) , T( 5,-5) , T( 6,-6) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { Y1.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    Y2.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T( 2,-2) , T(-1, 1) , T( 4,-4) , T(-3, 3) ,
               T( 6,-6) , T(-5, 5) , T( 8,-8) , T(-7, 7) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { Y2.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    Y2.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
   { Y2.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    Y2.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { T(-2, 2) , T( 1,-1) , T(-4, 4) , T( 3,-3) ,
               T(-6, 6) , T( 5,-5) , T(-8, 8) , T( 7,-7) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { Y2.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    Y0.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { T(  9, -9) , T( 10,-10) , T( 11,-11) , T( 12,-12) ,
                 T( 13,-13) , T( 14,-14) , T( 15,-15) , T( 16,-16) ,
                 T( -1,  1) , T( -2,  2) , T( -3,  3) , T( -4,  4) ,
                 T( -5,  5) , T( -6,  6) , T( -7,  7) , T( -8,  8) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { Y0.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    Y0.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
   { Y0.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    Y0.apply( qclab::Op::Trans , 4 , vec4 ) ;
    check4 = { T( -9,  9) , T(-10, 10) , T(-11, 11) , T(-12, 12) ,
               T(-13, 13) , T(-14, 14) , T(-15, 15) , T(-16, 16) ,
               T(  1, -1) , T(  2, -2) , T(  3, -3) , T(  4, -4) ,
               T(  5, -5) , T(  6, -6) , T(  7, -7) , T(  8, -8) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { Y0.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    Y1.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = { T(  5, -5) , T(  6, -6) , T(  7, -7) , T(  8, -8) ,
               T( -1,  1) , T( -2,  2) , T( -3,  3) , T( -4,  4) ,
               T( 13,-13) , T( 14,-14) , T( 15,-15) , T( 16,-16) ,
               T( -9,  9) , T(-10, 10) , T(-11, 11) , T(-12, 12) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { Y1.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    Y1.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
   { Y1.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    Y1.apply( qclab::Op::Trans , 4 , vec4 ) ;
    check4 = { T( -5,  5) , T( -6,  6) , T( -7,  7) , T( -8,  8) ,
               T(  1, -1) , T(  2, -2) , T(  3, -3) , T(  4, -4) ,
               T(-13, 13) , T(-14, 14) , T(-15, 15) , T(-16, 16) ,
               T(  9, -9) , T( 10,-10) , T( 11,-11) , T( 12,-12) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { Y1.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    Y2.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = { T(  3, -3) , T(  4, -4) , T( -1,  1) , T( -2,  2) ,
               T(  7, -7) , T(  8, -8) , T( -5,  5) , T( -6,  6) ,
               T( 11,-11) , T( 12,-12) , T( -9,  9) , T(-10, 10) ,
               T( 15,-15) , T( 16,-16) , T(-13, 13) , T(-14, 14) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { Y2.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    Y2.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
   { Y2.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    Y2.apply( qclab::Op::Trans , 4 , vec4 ) ;
    check4 = { T( -3,  3) , T( -4,  4) , T(  1, -1) , T(  2, -2) ,
               T( -7,  7) , T( -8,  8) , T(  5, -5) , T(  6, -6) ,
               T(-11, 11) , T(-12, 12) , T(  9, -9) , T( 10,-10) ,
               T(-15, 15) , T(-16, 16) , T( 13,-13) , T( 14,-14) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { Y2.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    Y3.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = { T(  2, -2) , T( -1,  1) , T(  4, -4) , T( -3,  3) ,
               T(  6, -6) , T( -5,  5) , T(  8, -8) , T( -7,  7) ,
               T( 10,-10) , T( -9,  9) , T( 12,-12) , T(-11, 11) ,
               T( 14,-14) , T(-13, 13) , T( 16,-16) , T(-15, 15) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { Y3.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    Y3.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
   { Y3.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    Y3.apply( qclab::Op::Trans , 4 , vec4 ) ;
    check4 = { T( -2,  2) , T(  1, -1) , T( -4,  4) , T(  3, -3) ,
               T( -6,  6) , T(  5, -5) , T( -8,  8) , T(  7, -7) ,
               T(-10, 10) , T(  9, -9) , T(-12, 12) , T( 11,-11) ,
               T(-14, 14) , T( 13,-13) , T(-16, 16) , T( 15,-15) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { Y3.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif
  }

}


/*
 * complex float
 */
TEST( qclab_qgates_PauliY , complex_float ) {
  test_qclab_qgates_PauliY< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_PauliY , complex_double ) {
  test_qclab_qgates_PauliY< std::complex< double > >() ;
}

