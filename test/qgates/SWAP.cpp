#include <gtest/gtest.h>
#include "qclab/qgates/SWAP.hpp"
#include "qclab/dense/kron.hpp"

template <typename T>
void test_qclab_qgates_SWAP() {

  const auto I1 = qclab::dense::eye< T >(  2 ) ;
  const auto I2 = qclab::dense::eye< T >(  4 ) ;
  const auto I3 = qclab::dense::eye< T >(  8 ) ;
  const auto I4 = qclab::dense::eye< T >( 16 ) ;

  using V = std::vector< T > ;
  const V v2 = { 3 , 5 , 2 , 7 } ;
  const V v3 = { 3 , 5 , 2 , 7 , 4 , 1 , 8 , 3 } ;
  const V v4 = { 3 , 5 , 2 , 7 , 4 , 1 , 8 , 3 , 7 , 2 , 5 , 6 , 8 , 9 , 5 , 1};

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
    auto vec2 = v2 ;
    swap.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 3 , 2 , 5 , 7 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { swap.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // apply (3 qubits)
    auto vec3 = v3 ;
    swap.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 3 , 5 , 4 , 1 , 2 , 7 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { swap.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    int qnew[] = { 1 , 2 } ;
    swap.setQubits( &qnew[0] ) ;
    vec3 = v3 ;
    swap.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , 2 , 5 , 7 , 4 , 8 , 1 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { swap.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // apply (4 qubits)
    auto vec4 = v4 ;
    swap.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { 3 , 5 , 4 , 1 , 2 , 7 , 8 , 3 , 7 , 2 , 8 , 9 , 5 , 6 , 5 , 1};
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { swap.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    qnew[0] = 0 ;
    qnew[1] = 1 ;
    swap.setQubits( &qnew[0] ) ;
    vec4 = v4 ;
    swap.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = { 3 , 5 , 2 , 7 , 7 , 2 , 5 , 6 , 4 , 1 , 8 , 3 , 8 , 9 , 5 , 1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { swap.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    qnew[0] = 2 ;
    qnew[1] = 3 ;
    swap.setQubits( &qnew[0] ) ;
    vec4 = v4 ;
    swap.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = { 3 , 2 , 5 , 7 , 4 , 8 , 1 , 3 , 7 , 5 , 2 , 6 , 8 , 5 , 9 , 1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { swap.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif
  }

  {
    qclab::qgates::SWAP< T >  swap( 0 , 2 ) ;

    // apply (3 qubits)
    auto vec3 = v3 ;
    swap.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 3 , 4 , 2 , 8 , 5 , 1 , 7 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { swap.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif
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

