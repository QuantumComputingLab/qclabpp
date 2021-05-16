#include <gtest/gtest.h>
#include "QRotation.hpp"

template <typename T>
void test_qclab_QRotation() {

  const T pi = 4 * std::atan(1) ;
  const T eps = std::numeric_limits< T >::epsilon() ;

  {
    qclab::QRotation< T >  R ;

    EXPECT_FALSE( R.fixed() ) ;         // fixed
    EXPECT_EQ( R.cos() , 1.0 ) ;        // cos
    EXPECT_EQ( R.sin() , 0.0 ) ;        // sin
    EXPECT_EQ( R.theta() , 0.0 ) ;      // theta

    // fixed
    R.makeFixed() ;
    EXPECT_TRUE( R.fixed() ) ;
    R.makeVariable() ;
    EXPECT_FALSE( R.fixed() ) ;

    // update(angle)
    qclab::QAngle< T >  new_angle( 0.5 ) ;
    R.update( new_angle ) ;
    EXPECT_NEAR( R.theta() , 1 , eps ) ;
    EXPECT_NEAR( R.cos() , std::cos( 0.5 ) , eps ) ;
    EXPECT_NEAR( R.sin() , std::sin( 0.5 ) , eps ) ;

    // update(theta)
    R.update( pi/2 ) ;
    EXPECT_NEAR( R.theta() , pi/2 , eps ) ;
    EXPECT_NEAR( R.cos() , std::cos( pi/4 ) , eps ) ;
    EXPECT_NEAR( R.sin() , std::sin( pi/4 ) , eps ) ;

    // update(cos,sin)
    R.update( std::cos( pi/3 ) , std::sin( pi/3 ) ) ;
    EXPECT_NEAR( R.cos() , std::cos( pi/3 ) , eps ) ;
    EXPECT_NEAR( R.sin() , std::sin( pi/3 ) , eps ) ;
    EXPECT_NEAR( R.theta() , 2*(pi/3) , eps ) ;

    // operators == and !=
    qclab::QRotation< T >  R2( std::cos( pi/3 ) , std::sin( pi/3 ) ) ;
    EXPECT_TRUE( R == R2 ) ;
    EXPECT_FALSE( R != R2 ) ;
    R2.update( 1 ) ;
    EXPECT_TRUE( R != R2 ) ;
    EXPECT_FALSE( R == R2 ) ;
}

  {
    qclab::QAngle< T > angle( pi/2 ) ;
    qclab::QRotation< T >  R( angle ) ;

    EXPECT_FALSE( R.fixed() ) ;                        // fixed
    EXPECT_NEAR( R.theta() , pi , eps ) ;              // theta
    EXPECT_NEAR( R.cos() , std::cos( pi/2 ) , eps ) ;  // cos
    EXPECT_NEAR( R.sin() , std::sin( pi/2 ) , eps ) ;  // sin
  }

  {
    qclab::QAngle< T > angle( pi/2 ) ;
    qclab::QRotation< T >  R( angle , true ) ;

    EXPECT_TRUE( R.fixed() ) ;                         // fixed
    EXPECT_NEAR( R.theta() , pi , eps ) ;              // theta
    EXPECT_NEAR( R.cos() , std::cos( pi/2 ) , eps ) ;  // cos
    EXPECT_NEAR( R.sin() , std::sin( pi/2 ) , eps ) ;  // sin
  }

  {
    qclab::QRotation< T >  R( pi/2 ) ;

    EXPECT_FALSE( R.fixed() ) ;                        // fixed
    EXPECT_NEAR( R.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( R.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( R.sin() , std::sin( pi/4 ) , eps ) ;  // sin
  }

  {
    qclab::QRotation< T >  R( pi/2 , true ) ;

    EXPECT_TRUE( R.fixed() ) ;                         // fixed
    EXPECT_NEAR( R.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( R.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( R.sin() , std::sin( pi/4 ) , eps ) ;  // sin
  }

  {
    const T cos = std::cos( pi/4 ) ;
    const T sin = std::sin( pi/4 ) ;
    qclab::QRotation< T >  R( cos , sin ) ;

    EXPECT_FALSE( R.fixed() ) ;              // fixed
    EXPECT_NEAR( R.theta() , pi/2 , eps ) ;  // theta
    EXPECT_NEAR( R.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( R.sin() , sin , eps ) ;     // sin
  }

  {
    const T cos = std::cos( pi/4 ) ;
    const T sin = std::sin( pi/4 ) ;
    qclab::QRotation< T >  R( cos , sin , true ) ;

    EXPECT_TRUE( R.fixed() ) ;               // fixed
    EXPECT_NEAR( R.theta() , pi/2 , eps ) ;  // theta
    EXPECT_NEAR( R.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( R.sin() , sin , eps ) ;     // sin
  }

}


/*
 * float
 */
TEST( qclab_QRotation , float ) {
  test_qclab_QRotation< float >() ;
}

/*
 * double
 */
TEST( qclab_QRotation , double ) {
  test_qclab_QRotation< double >() ;
}

