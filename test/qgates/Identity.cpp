#include <gtest/gtest.h>
#include "qclab/qgates/Identity.hpp"

template <typename T>
void test_qclab_qgates_Identity() {

  qclab::qgates::Identity< T >  I ;

  EXPECT_EQ( I.nbQubits() , 1 ) ;   // nbQubits
  EXPECT_TRUE( I.fixed() ) ;        // fixed
  EXPECT_FALSE( I.controlled() ) ;  // controlled

  // qubit
  EXPECT_EQ( I.qubit() , 0 ) ;
  I.setQubit( 2 ) ;
  EXPECT_EQ( I.qubit() , 2 ) ;

  // qubits
  auto qubits = I.qubits() ;
  EXPECT_EQ( qubits.size() , 1 ) ;
  EXPECT_EQ( qubits[0] , 2 ) ;
  int qnew = 3 ;
  I.setQubits( &qnew ) ;
  EXPECT_EQ( I.qubit() , 3 ) ;

  // matrix
  EXPECT_EQ( I.matrix()(0,0) , T(1.0) ) ;
  EXPECT_EQ( I.matrix()(1,0) , T(0.0) ) ;
  EXPECT_EQ( I.matrix()(0,1) , T(0.0) ) ;
  EXPECT_EQ( I.matrix()(1,1) , T(1.0) ) ;

  // print
  I.print() ;

  // toQASM
  std::stringstream qasm ;
  EXPECT_EQ( I.toQASM( qasm ) , 0 ) ;
  EXPECT_EQ( qasm.str() , "" ) ;

  // operators == and !=
  qclab::qgates::Identity< T > I2 ;
  EXPECT_TRUE( I == I2 ) ;
  EXPECT_FALSE( I != I2 ) ;

  // apply
  {
    using V = std::vector< T > ;
    const V v1 = { 0 , 1 } ;

    qclab::qgates::Identity< T >  I0( 0 ) ;

    // nbQubits = 1
    auto vec1 = v1 ;
    I0.apply( qclab::Op::NoTrans , 1 , vec1 ) ;
    V check1 = { 0 , 1 } ;
    EXPECT_TRUE( vec1 == check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ; T* vec1_ = vec1.data() ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { I0.apply_device( qclab::Op::NoTrans , 1 , vec1_ ) ; }
    EXPECT_TRUE( vec1 == check1 ) ;
  #endif

    vec1 = v1 ;
    I0.apply( qclab::Op::Trans , 1 , vec1 ) ;
    EXPECT_TRUE( vec1 == check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { I0.apply_device( qclab::Op::Trans , 1 , vec1_ ) ; }
    EXPECT_TRUE( vec1 == check1 ) ;
  #endif

    vec1 = v1 ;
    I0.apply( qclab::Op::ConjTrans , 1 , vec1 ) ;
    EXPECT_TRUE( vec1 == check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { I0.apply_device( qclab::Op::ConjTrans , 1 , vec1_ ) ; }
    EXPECT_TRUE( vec1 == check1 ) ;
  #endif
  }

}


/*
 * float
 */
TEST( qclab_qgates_Identity , float ) {
  test_qclab_qgates_Identity< float >() ;
}

/*
 * double
 */
TEST( qclab_qgates_Identity , double ) {
  test_qclab_qgates_Identity< double >() ;
}


/*
 * complex float
 */
TEST( qclab_qgates_Identity , complex_float ) {
  test_qclab_qgates_Identity< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_Identity , complex_double ) {
  test_qclab_qgates_Identity< std::complex< double > >() ;
}

