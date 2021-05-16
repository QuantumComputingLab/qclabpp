#include <gtest/gtest.h>
#include "dense/transpose.hpp"

template <typename T>
void test_qclab_dense_transpose() {

  using M = qclab::dense::SquareMatrix< T > ;

  M matrix1( 1 , 2 ,
             3 , 4 ) ;
  M matrix2( 1 , 5 ,  9 , 13 ,
             2 , 6 , 10 , 14 ,
             3 , 7 , 11 , 15 ,
             4 , 8 , 12 , 16 ) ;

  M check1( 1 , 2 ,
            3 , 4 ) ;
  M check2( 1 , 5 ,  9 , 13 ,
            2 , 6 , 10 , 14 ,
            3 , 7 , 11 , 15 ,
            4 , 8 , 12 , 16 ) ;

  M check1T( 1 , 3 ,
             2 , 4 ) ;
  M check2T(  1 ,  2 ,  3 ,  4 ,
              5 ,  6 ,  7 ,  8 ,
              9 , 10 , 11 , 12 ,
             13 , 14 , 15 , 16 ) ;

  // NoTrans
  {
    M m1( matrix1 ) ;
    qclab::dense::operateInPlace( qclab::Op::NoTrans , m1 ) ;
    EXPECT_EQ( m1 , check1 ) ;
    M m2( matrix2 ) ;
    qclab::dense::operateInPlace( qclab::Op::NoTrans , m2 ) ;
    EXPECT_EQ( m2 , check2 ) ;
    M m3 = qclab::dense::operate( qclab::Op::NoTrans , matrix2 ) ;
    EXPECT_EQ( m3 , check2 ) ;
  }

  // Trans
  {
    M m1( matrix1 ) ;
    qclab::dense::operateInPlace( qclab::Op::Trans , m1 ) ;
    EXPECT_EQ( m1 , check1T ) ;
    M m2( matrix2 ) ;
    qclab::dense::operateInPlace( qclab::Op::Trans , m2 ) ;
    EXPECT_EQ( m2 , check2T ) ;
    M m3 = qclab::dense::operate( qclab::Op::Trans , matrix2 ) ;
    EXPECT_EQ( m3 , check2T ) ;
    M m4( matrix2 ) ;
    qclab::dense::transInPlace( m4 ) ;
    EXPECT_EQ( m4 , check2T ) ;
    M m5 = qclab::dense::trans( matrix2 ) ;
    EXPECT_EQ( m5 , check2T ) ;
  }

  // ConjTrans
  {
    M m1( matrix1 ) ;
    qclab::dense::operateInPlace( qclab::Op::ConjTrans , m1 ) ;
    EXPECT_EQ( m1 , check1T ) ;
    M m2( matrix2 ) ;
    qclab::dense::operateInPlace( qclab::Op::ConjTrans , m2 ) ;
    EXPECT_EQ( m2 , check2T ) ;
    M m3 = qclab::dense::operate( qclab::Op::ConjTrans , matrix2 ) ;
    EXPECT_EQ( m3 , check2T ) ;
    M m4( matrix2 ) ;
    qclab::dense::conjTransInPlace( m4 ) ;
    EXPECT_EQ( m4 , check2T ) ;
    M m5 = qclab::dense::conjTrans( matrix2 ) ;
    EXPECT_EQ( m5 , check2T ) ;
  }

}


template <typename T>
void test_qclab_dense_transpose_complex() {

  using M = qclab::dense::SquareMatrix< T > ;

  M matrix1( T(1,5) , T(2,6) ,
             T(3,7) , T(4,8) ) ;
  M matrix2( T(0,1) , T(0,5) , T(0, 9) , T(0,13) ,
             T(0,2) , T(0,6) , T(0,10) , T(0,14) ,
             T(0,3) , T(0,7) , T(0,11) , T(0,15) ,
             T(0,4) , T(0,8) , T(0,12) , T(0,16) ) ;

  M check1( T(1,5) , T(2,6) ,
            T(3,7) , T(4,8) ) ;
  M check2( T(0,1) , T(0,5) , T(0, 9) , T(0,13) ,
            T(0,2) , T(0,6) , T(0,10) , T(0,14) ,
            T(0,3) , T(0,7) , T(0,11) , T(0,15) ,
            T(0,4) , T(0,8) , T(0,12) , T(0,16) ) ;

  M check1T( T(1,5) , T(3,7) ,
             T(2,6) , T(4,8) ) ;
  M check2T( T(0, 1) , T(0, 2) , T(0, 3) , T(0, 4) ,
             T(0, 5) , T(0, 6) , T(0, 7) , T(0, 8) ,
             T(0, 9) , T(0,10) , T(0,11) , T(0,12) ,
             T(0,13) , T(0,14) , T(0,15) , T(0,16) ) ;

  M check1H( T(1,-5) , T(3,-7) ,
             T(2,-6) , T(4,-8) ) ;
  M check2H( T(0, -1) , T(0, -2) , T(0, -3) , T(0, -4) ,
             T(0, -5) , T(0, -6) , T(0, -7) , T(0, -8) ,
             T(0, -9) , T(0,-10) , T(0,-11) , T(0,-12) ,
             T(0,-13) , T(0,-14) , T(0,-15) , T(0,-16) ) ;

  // NoTrans
  {
    M m1( matrix1 ) ;
    qclab::dense::operateInPlace( qclab::Op::NoTrans , m1 ) ;
    EXPECT_EQ( m1 , check1 ) ;
    M m2( matrix2 ) ;
    qclab::dense::operateInPlace( qclab::Op::NoTrans , m2 ) ;
    EXPECT_EQ( m2 , check2 ) ;
    M m3 = qclab::dense::operate( qclab::Op::NoTrans , matrix2 ) ;
    EXPECT_EQ( m3 , check2 ) ;
  }

  // Trans
  {
    M m1( matrix1 ) ;
    qclab::dense::operateInPlace( qclab::Op::Trans , m1 ) ;
    EXPECT_EQ( m1 , check1T ) ;
    M m2( matrix2 ) ;
    qclab::dense::operateInPlace( qclab::Op::Trans , m2 ) ;
    EXPECT_EQ( m2 , check2T ) ;
    M m3 = qclab::dense::operate( qclab::Op::Trans , matrix2 ) ;
    EXPECT_EQ( m3 , check2T ) ;
    M m4( matrix2 ) ;
    qclab::dense::transInPlace( m4 ) ;
    EXPECT_EQ( m4 , check2T ) ;
    M m5 = qclab::dense::trans( matrix2 ) ;
    EXPECT_EQ( m5 , check2T ) ;
  }

  // ConjTrans
  {
    M m1( matrix1 ) ;
    qclab::dense::operateInPlace( qclab::Op::ConjTrans , m1 ) ;
    EXPECT_EQ( m1 , check1H ) ;
    M m2( matrix2 ) ;
    qclab::dense::operateInPlace( qclab::Op::ConjTrans , m2 ) ;
    EXPECT_EQ( m2 , check2H ) ;
    M m3 = qclab::dense::operate( qclab::Op::ConjTrans , matrix2 ) ;
    EXPECT_EQ( m3 , check2H ) ;
    M m4( matrix2 ) ;
    qclab::dense::conjTransInPlace( m4 ) ;
    EXPECT_EQ( m4 , check2H ) ;
    M m5 = qclab::dense::conjTrans( matrix2 ) ;
    EXPECT_EQ( m5 , check2H ) ;
  }

}


/*
 * float
 */
TEST( qclab_dense_transpose , float ) {
  test_qclab_dense_transpose< float >() ;
}

/*
 * double
 */
TEST( qclab_dense_transpose , double ) {
  test_qclab_dense_transpose< double >() ;
}


/*
 * complex float
 */
TEST( qclab_dense_transpose , complex_float ) {
  test_qclab_dense_transpose< std::complex< float > >() ;
  test_qclab_dense_transpose_complex< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_dense_transpose , complex_double ) {
  test_qclab_dense_transpose< std::complex< double > >() ;
  test_qclab_dense_transpose_complex< std::complex< double > >() ;
}

