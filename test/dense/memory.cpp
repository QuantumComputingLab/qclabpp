#include <gtest/gtest.h>
#include "dense/memory.hpp"
#include <complex>

template <typename T>
void test_qclab_dense_memory() {

  // alloc_unique_array
  std::unique_ptr< T[] >  alloc1 = qclab::dense::alloc_unique_array< T >( 1 ) ;
  EXPECT_EQ( alloc1[0] , alloc1[0] ) ;

  std::unique_ptr< T[] >  alloc2 = qclab::dense::alloc_unique_array< T >( 2 ) ;
  EXPECT_EQ( alloc2[0] , alloc2[0] ) ;
  EXPECT_EQ( alloc2[1] , alloc2[1] ) ;

  std::unique_ptr< T[] >  alloc10 = qclab::dense::alloc_unique_array< T >( 10 );
  for ( int64_t i = 0; i < 1; i++ ) {
    EXPECT_EQ( alloc10[i] , alloc10[i] ) ;
  }

  // init_unique_array (zero)
  std::unique_ptr< T[] >  mem1 = qclab::dense::init_unique_array< T >( 1 ) ;
  EXPECT_EQ( mem1[0] , T(0) ) ;

  std::unique_ptr< T[] >  mem2 = qclab::dense::init_unique_array< T >( 2 ) ;
  EXPECT_EQ( mem2[0] , T(0) ) ;
  EXPECT_EQ( mem2[1] , T(0) ) ;

  std::unique_ptr< T[] >  mem10 = qclab::dense::init_unique_array< T >( 10 ) ;
  for ( int64_t i = 0; i < 1; i++ ) {
    EXPECT_EQ( mem10[i] , T(0) ) ;
  }

  // init_unique_array (zero)
  auto init1 = qclab::dense::init_unique_array< T >( 1 , T(3.14) ) ;
  EXPECT_EQ( init1[0] , T(3.14) ) ;

  auto init2 = qclab::dense::init_unique_array< T >( 2 , T(3.14) ) ;
  EXPECT_EQ( init2[0] , T(3.14) ) ;
  EXPECT_EQ( init2[1] , T(3.14) ) ;

  auto init10 = qclab::dense::init_unique_array< T >( 10 , T(3.14) ) ;
  for ( int64_t i = 0; i < 1; i++ ) {
    EXPECT_EQ( init10[i] , T(3.14) ) ;
  }

}


/*
 * float
 */
TEST( qclab_dense_memory , float ) {
  test_qclab_dense_memory< float >() ;
}

/*
 * double
 */
TEST( qclab_dense_memory , double ) {
  test_qclab_dense_memory< double >() ;
}


/*
 * complex float
 */
TEST( qclab_dense_memory , complex_float ) {
  test_qclab_dense_memory< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_dense_memory , complex_double ) {
  test_qclab_dense_memory< std::complex< double > >() ;
}

