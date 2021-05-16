#include <gtest/gtest.h>
#include "dense/SquareMatrix.hpp"
#include <complex>

template <typename T>
void test_qclab_dense_SquareMatrix() {

  qclab::dense::SquareMatrix< T >  M1( 3 ) ;
  EXPECT_EQ( M1.size() , 3 ) ;
  EXPECT_EQ( M1.rows() , 3 ) ;
  EXPECT_EQ( M1.cols() , 3 ) ;
  EXPECT_EQ( M1.ld() , 3 ) ;

  int c = 0 ;
  for ( int i = 0; i < 3; i++ ) {
    for ( int j = 0; j < 3; j++ ) {
      M1(i,j) = c + 0.5 ;
      c++ ;
    }
  }
  c = 0 ;
  for ( int i = 0; i < 3; i++ ) {
    for ( int j = 0; j < 3; j++ ) {
      EXPECT_EQ( M1(i,j) , T(c + 0.5) ) ;
      c++ ;
    }
  }

  qclab::dense::SquareMatrix< T >  M2( 3 , 3.14 ) ;
  EXPECT_EQ( M2.size() , 3 ) ;
  EXPECT_EQ( M2.rows() , 3 ) ;
  EXPECT_EQ( M2.cols() , 3 ) ;
  EXPECT_EQ( M2.ld() , 3 ) ;

  for ( int i = 0; i < 3; i++ ) {
    for ( int j = 0; j < 3; j++ ) {
      EXPECT_EQ( M2(i,j) , T(3.14) ) ;
    }
  }

  qclab::dense::SquareMatrix< T >  M3( 1.0 , 2.0 ,
                                       3.0 , 4.0 ) ;
  EXPECT_EQ( M3.size() , 2 ) ;
  EXPECT_EQ( M3.rows() , 2 ) ;
  EXPECT_EQ( M3.cols() , 2 ) ;
  EXPECT_EQ( M3.ld() , 2 ) ;

  EXPECT_EQ( M3(0,0) , T(1.0) ) ;
  EXPECT_EQ( M3(1,0) , T(3.0) ) ;
  EXPECT_EQ( M3(0,1) , T(2.0) ) ;
  EXPECT_EQ( M3(1,1) , T(4.0) ) ;

  qclab::dense::SquareMatrix< T >  M4(  1.0 ,  2.0 ,  3.0 ,  4.0 ,
                                        5.0 ,  6.0 ,  7.0 ,  8.0 ,
                                        9.0 , 10.0 , 11.0 , 12.0 ,
                                       13.0 , 14.0 , 15.0 , 16.0 ) ;
  EXPECT_EQ( M4.size() , 4 ) ;
  EXPECT_EQ( M4.rows() , 4 ) ;
  EXPECT_EQ( M4.cols() , 4 ) ;
  EXPECT_EQ( M4.ld() , 4 ) ;

  // column 1
  EXPECT_EQ( M4(0,0) , T( 1.0) ) ;
  EXPECT_EQ( M4(1,0) , T( 5.0) ) ;
  EXPECT_EQ( M4(2,0) , T( 9.0) ) ;
  EXPECT_EQ( M4(3,0) , T(13.0) ) ;
  // column 2
  EXPECT_EQ( M4(0,1) , T( 2.0) ) ;
  EXPECT_EQ( M4(1,1) , T( 6.0) ) ;
  EXPECT_EQ( M4(2,1) , T(10.0) ) ;
  EXPECT_EQ( M4(3,1) , T(14.0) ) ;
  // column 3
  EXPECT_EQ( M4(0,2) , T( 3.0) ) ;
  EXPECT_EQ( M4(1,2) , T( 7.0) ) ;
  EXPECT_EQ( M4(2,2) , T(11.0) ) ;
  EXPECT_EQ( M4(3,2) , T(15.0) ) ;
  // column 4
  EXPECT_EQ( M4(0,3) , T( 4.0) ) ;
  EXPECT_EQ( M4(1,3) , T( 8.0) ) ;
  EXPECT_EQ( M4(2,3) , T(12.0) ) ;
  EXPECT_EQ( M4(3,3) , T(16.0) ) ;

  qclab::dense::SquareMatrix< T >  M4copy( M4 ) ;
  EXPECT_EQ( M4copy.size() , 4 ) ;
  EXPECT_EQ( M4copy.rows() , 4 ) ;
  EXPECT_EQ( M4copy.cols() , 4 ) ;
  EXPECT_EQ( M4copy.ld() , 4 ) ;
  EXPECT_TRUE( M4copy == M4 ) ;

  // operator =
  qclab::dense::SquareMatrix< T >  M1copy ;
  M1copy = M1 ;
  EXPECT_EQ( M1copy.size() , 3 ) ;
  EXPECT_EQ( M1copy.rows() , 3 ) ;
  EXPECT_EQ( M1copy.cols() , 3 ) ;
  EXPECT_EQ( M1copy.ld() , 3 ) ;
  EXPECT_TRUE( M1copy == M1 ) ;

  qclab::dense::SquareMatrix< T >  M2copy( 3 ) ;
  M2copy = M2 ;
  EXPECT_EQ( M2copy.size() , 3 ) ;
  EXPECT_EQ( M2copy.rows() , 3 ) ;
  EXPECT_EQ( M2copy.cols() , 3 ) ;
  EXPECT_EQ( M2copy.ld() , 3 ) ;
  EXPECT_TRUE( M2copy == M2 ) ;

  // operators == and !=
  EXPECT_TRUE( M2 == M2 ) ; EXPECT_FALSE( M2 != M2 ) ;
  EXPECT_TRUE( M3 == M3 ) ; EXPECT_FALSE( M3 != M3 ) ;
  EXPECT_TRUE( M4 == M4 ) ; EXPECT_FALSE( M4 != M4 ) ;
  EXPECT_TRUE( M2 != M3 ) ; EXPECT_FALSE( M2 == M3 ) ;
  EXPECT_TRUE( M3 != M4 ) ; EXPECT_FALSE( M3 == M4 ) ;

  qclab::dense::SquareMatrix< T >  M5(  1.0 ,  2.0 ,  3.0 ,  4.0 ,
                                        5.0 ,  6.0 ,  0.0 ,  8.0 ,
                                        9.0 , 10.0 , 11.0 , 12.0 ,
                                       13.0 , 14.0 , 15.0 , 16.0 ) ;
  EXPECT_TRUE( M4 != M5 ) ; EXPECT_FALSE( M4 == M5 ) ;
  M5(1,2) = 7.0 ;
  EXPECT_TRUE( M4 == M5 ) ; EXPECT_FALSE( M4 != M5 ) ;

  // operators +=, -=, *=, +, - and *
  qclab::dense::SquareMatrix< T >  A( 1 , 2 ,
                                      3 , 4 ) ;
  qclab::dense::SquareMatrix< T >  B( 5 , 7 ,
                                      6 , 8 ) ;
  qclab::dense::SquareMatrix< T >  C(  9 , 11 ,
                                      10 , 12 ) ;
  qclab::dense::SquareMatrix< T >  D( 13 , 14 ,
                                      15 , 16 ) ;

  qclab::dense::SquareMatrix< T >  sumAB_check(  6 ,  9 ,
                                                 9 , 12 ) ;
  qclab::dense::SquareMatrix< T >  sumABC_check( 15 , 20 ,
                                                 19 , 24 ) ;
  qclab::dense::SquareMatrix< T >  sumABCD_check( 28 , 34 ,
                                                  34 , 40 ) ;

  auto sumAB = A + B ;
  EXPECT_EQ( sumAB.size() , 2 ) ;
  EXPECT_TRUE( sumAB == sumAB_check ) ;

  auto sumABC = A + B + C ;
  EXPECT_EQ( sumABC.size() , 2 ) ;
  EXPECT_TRUE( sumABC == sumABC_check ) ;

  auto sumABCD = A + B + C + D ;
  EXPECT_EQ( sumABCD.size() , 2 ) ;
  EXPECT_TRUE( sumABCD == sumABCD_check ) ;

  auto sum = A ;
  sum += B ;  EXPECT_TRUE( sum == sumAB_check ) ;
  sum += C ;  EXPECT_TRUE( sum == sumABC_check ) ;
  sum += D ;  EXPECT_TRUE( sum == sumABCD_check ) ;

  qclab::dense::SquareMatrix< T >  minusAB_check( -4 , -5 ,
                                                  -3 , -4 ) ;
  qclab::dense::SquareMatrix< T >  minusABC_check( -13 , -16 ,
                                                   -13 , -16 ) ;
  qclab::dense::SquareMatrix< T >  minusABCD_check( -26 , -30 ,
                                                    -28 , -32 ) ;

  auto minusAB = A - B ;
  EXPECT_EQ( minusAB.size() , 2 ) ;
  EXPECT_TRUE( minusAB == minusAB_check ) ;

  auto minusABC = A - B - C ;
  EXPECT_EQ( minusABC.size() , 2 ) ;
  EXPECT_TRUE( minusABC == minusABC_check ) ;

  auto minusABCD = A - B - C - D ;
  EXPECT_EQ( minusABCD.size() , 2 ) ;
  EXPECT_TRUE( minusABCD == minusABCD_check ) ;

  auto minus = A ;
  minus -= B ;  EXPECT_TRUE( minus == minusAB_check ) ;
  minus -= C ;  EXPECT_TRUE( minus == minusABC_check ) ;
  minus -= D ;  EXPECT_TRUE( minus == minusABCD_check ) ;

  qclab::dense::SquareMatrix< T >  prodAB_check( 17 , 23 ,
                                                 39 , 53 ) ;
  qclab::dense::SquareMatrix< T >  prodABC_check( 383 ,  463 ,
                                                  881 , 1065 ) ;
  qclab::dense::SquareMatrix< T >  prodABCD_check( 11924 , 12770 ,
                                                   27428 , 29374 ) ;

  auto prodAB = A * B ;
  EXPECT_EQ( prodAB.size() , 2 ) ;
  EXPECT_TRUE( prodAB == prodAB_check ) ;

  auto prodABC = A * B * C ;
  EXPECT_EQ( prodABC.size() , 2 ) ;
  EXPECT_TRUE( prodABC == prodABC_check ) ;

  auto prodABCD = A * B * C * D ;
  EXPECT_EQ( prodABCD.size() , 2 ) ;
  EXPECT_TRUE( prodABCD == prodABCD_check ) ;

  auto prod = A ;
  prod *= B ;  EXPECT_TRUE( prod == prodAB_check ) ;
  prod *= C ;  EXPECT_TRUE( prod == prodABC_check ) ;
  prod *= D ;  EXPECT_TRUE( prod == prodABCD_check ) ;

  // zeros
  auto Z = qclab::dense::zeros< T >( 2 ) ;
  qclab::dense::SquareMatrix< T >  Z_check( 0 , 0 ,
                                            0 , 0 ) ;
  EXPECT_TRUE( Z == Z_check ) ;

  // eye
  auto I = qclab::dense::eye< T >( 2 ) ;
  qclab::dense::SquareMatrix< T >  I_check( 1 , 0 ,
                                            0 , 1 ) ;
  EXPECT_TRUE( I == I_check ) ;

}


/*
 * float
 */
TEST( qclab_dense_SquareMatrix , float ) {
  test_qclab_dense_SquareMatrix< float >() ;
}

/*
 * double
 */
TEST( qclab_dense_SquareMatrix , double ) {
  test_qclab_dense_SquareMatrix< double >() ;
}


/*
 * complex float
 */
TEST( qclab_dense_SquareMatrix , complex_float ) {
  test_qclab_dense_SquareMatrix< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_dense_SquareMatrix , complex_double ) {
  test_qclab_dense_SquareMatrix< std::complex< double > >() ;
}

