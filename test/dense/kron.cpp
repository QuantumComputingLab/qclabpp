#include <gtest/gtest.h>
#include "qclab/dense/kron.hpp"
#include <complex>

template <typename T>
void test_qclab_dense_kron() {

  qclab::dense::SquareMatrix< T >  I( 1 , 0 ,
                                      0 , 1 ) ;
  qclab::dense::SquareMatrix< T >  X( 0 , 1 ,
                                      1 , 0 ) ;
//std::cout << "I =\n" ; qclab::printMatrix2x2( I ) ;
//std::cout << "X =\n" ; qclab::printMatrix2x2( X ) ;

  qclab::dense::SquareMatrix< T >  IoX = qclab::dense::kron( I , X ) ;
  qclab::dense::SquareMatrix< T >  IoX_check( 0 , 1 , 0 , 0 ,
                                              1 , 0 , 0 , 0 ,
                                              0 , 0 , 0 , 1 ,
                                              0 , 0 , 1 , 0 ) ;
  EXPECT_TRUE( IoX == IoX_check ) ;
//std::cout << "I otimes X =\n" ; qclab::printMatrix4x4( IoX ) ;

  qclab::dense::SquareMatrix< T >  XoI = qclab::dense::kron( X , I ) ;
  qclab::dense::SquareMatrix< T >  XoI_check( 0 , 0 , 1 , 0 ,
                                              0 , 0 , 0 , 1 ,
                                              1 , 0 , 0 , 0 ,
                                              0 , 1 , 0 , 0 ) ;
  EXPECT_TRUE( XoI == XoI_check ) ;
//std::cout << "X otimes I =\n" ; qclab::printMatrix4x4( XoI ) ;

  qclab::dense::SquareMatrix< T >  A(  1 ,  2 ,
                                       3 ,  5 ) ;
  qclab::dense::SquareMatrix< T >  B(  7 , 11 ,
                                      13 , 17 ) ;
//std::cout << "A =\n" ; qclab::printMatrix2x2( A ) ;
//std::cout << "B =\n" ; qclab::printMatrix2x2( B ) ;

  qclab::dense::SquareMatrix< T >  AoB = qclab::dense::kron( A , B ) ;
  qclab::dense::SquareMatrix< T >  AoB_check(  7 , 11 , 14 , 22 ,
                                              13 , 17 , 26 , 34 ,
                                              21 , 33 , 35 , 55 ,
                                              39 , 51 , 65 , 85 ) ;
  EXPECT_TRUE( AoB == AoB_check ) ;
//std::cout << "A otimes B =\n" ; qclab::printMatrix4x4( AoB ) ;

  qclab::dense::SquareMatrix< T >  BoA = qclab::dense::kron( B , A ) ;
  qclab::dense::SquareMatrix< T >  BoA_check(  7 , 14 , 11 , 22 ,
                                              21 , 35 , 33 , 55 ,
                                              13 , 26 , 17 , 34 ,
                                              39 , 65 , 51 , 85 ) ;
  EXPECT_TRUE( BoA == BoA_check ) ;
//std::cout << "B otimes A =\n" ; qclab::printMatrix4x4( BoA ) ;

}


/*
 * float
 */
TEST( qclab_dense_kron , float ) {
  test_qclab_dense_kron< float >() ;
}

/*
 * double
 */
TEST( qclab_dense_kron , double ) {
  test_qclab_dense_kron< double >() ;
}


/*
 * complex float
 */
TEST( qclab_dense_kron , complex_float ) {
  test_qclab_dense_kron< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_dense_kron , complex_double ) {
  test_qclab_dense_kron< std::complex< double > >() ;
}

