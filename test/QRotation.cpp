#include <gtest/gtest.h>
#include "qclab/QRotation.hpp"

template <typename T>
void test_qclab_QRotation() {

  const T pi = 4 * std::atan(1) ;
  const T eps = std::numeric_limits< T >::epsilon() ;

  {
    qclab::QRotation< T >  R ;

    EXPECT_EQ( R.cos() , 1.0 ) ;    // cos
    EXPECT_EQ( R.sin() , 0.0 ) ;    // sin
    EXPECT_EQ( R.theta() , 0.0 ) ;  // theta

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

    EXPECT_NEAR( R.theta() , pi , eps ) ;              // theta
    EXPECT_NEAR( R.cos() , std::cos( pi/2 ) , eps ) ;  // cos
    EXPECT_NEAR( R.sin() , std::sin( pi/2 ) , eps ) ;  // sin
  }

  {
    qclab::QAngle< T > angle( pi/2 ) ;
    qclab::QRotation< T >  R( angle ) ;

    EXPECT_NEAR( R.theta() , pi , eps ) ;              // theta
    EXPECT_NEAR( R.cos() , std::cos( pi/2 ) , eps ) ;  // cos
    EXPECT_NEAR( R.sin() , std::sin( pi/2 ) , eps ) ;  // sin
  }

  {
    qclab::QRotation< T >  R( pi/2 ) ;

    EXPECT_NEAR( R.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( R.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( R.sin() , std::sin( pi/4 ) , eps ) ;  // sin
  }

  {
    qclab::QRotation< T >  R( pi/2 ) ;

    EXPECT_NEAR( R.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( R.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( R.sin() , std::sin( pi/4 ) , eps ) ;  // sin
  }

  {
    const T cos = std::cos( pi/4 ) ;
    const T sin = std::sin( pi/4 ) ;
    qclab::QRotation< T >  R( cos , sin ) ;

    EXPECT_NEAR( R.theta() , pi/2 , eps ) ;  // theta
    EXPECT_NEAR( R.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( R.sin() , sin , eps ) ;     // sin
  }

  {
    const T cos = std::cos( pi/4 ) ;
    const T sin = std::sin( pi/4 ) ;
    qclab::QRotation< T >  R( cos , sin ) ;

    EXPECT_NEAR( R.theta() , pi/2 , eps ) ;  // theta
    EXPECT_NEAR( R.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( R.sin() , sin , eps ) ;     // sin
  }

  {
    const T theta1 = pi/2 ;
    const qclab::QAngle< T >  angle1( theta1 ) ;
    qclab::QRotation< T >  R1( angle1 ) ;

    const T theta2 = pi/3 ;
    const qclab::QAngle< T >  angle2( theta2 ) ;
    qclab::QRotation< T >  R2( angle2 ) ;

    // operator *=
    R1 *= R2 ;
    qclab::QAngle< T >  angle = angle1 + angle2 ;
    T cos   = angle.cos() ;
    T sin   = angle.sin() ;
    T theta = angle.theta() ;
    EXPECT_NEAR( R1.cos() , cos , 10*eps ) ;                 // cos
    EXPECT_NEAR( R1.sin() , sin , 10*eps ) ;                 // sin
    EXPECT_NEAR( R1.theta() , 2*theta , 10*eps ) ;           // theta
    EXPECT_NEAR( R2.cos() , std::cos( theta2 ) , 10*eps ) ;  // cos
    EXPECT_NEAR( R2.sin() , std::sin( theta2 ) , 10*eps ) ;  // sin
    EXPECT_NEAR( R2.theta() , 2*theta2 , 10*eps ) ;          // theta
  }

  {
    const T theta1 = pi/2 ;
    const qclab::QAngle< T >  angle1( theta1 ) ;
    qclab::QRotation< T >  R1( angle1 ) ;

    const T theta2 = pi/3 ;
    const qclab::QAngle< T >  angle2( theta2 ) ;
    qclab::QRotation< T >  R2( angle2 ) ;

    // operator /=
    R1 /= R2 ;
    qclab::QAngle< T >  angle = angle1 - angle2 ;
    T cos   = angle.cos() ;
    T sin   = angle.sin() ;
    T theta = angle.theta() ;
    EXPECT_NEAR( R1.cos() , cos , 10*eps ) ;                 // cos
    EXPECT_NEAR( R1.sin() , sin , 10*eps ) ;                 // sin
    EXPECT_NEAR( R1.theta() , 2*theta , 10*eps ) ;           // theta
    EXPECT_NEAR( R2.cos() , std::cos( theta2 ) , 10*eps ) ;  // cos
    EXPECT_NEAR( R2.sin() , std::sin( theta2 ) , 10*eps ) ;  // sin
    EXPECT_NEAR( R2.theta() , 2*theta2 , 10*eps ) ;          // theta
  }

  {
    const T theta1 = pi/2 ;
    const qclab::QAngle< T >  angle1( theta1 ) ;
    qclab::QRotation< T >  R1( angle1 ) ;

    const T theta2 = pi/3 ;
    const qclab::QAngle< T >  angle2( theta2 ) ;
    qclab::QRotation< T >  R2( angle2 ) ;

    // operator *
    qclab::QRotation< T > R = R1 * R2 ;
    qclab::QAngle< T >  angle = angle1 + angle2 ;
    T cos   = angle.cos() ;
    T sin   = angle.sin() ;
    T theta = angle.theta() ;
    EXPECT_NEAR( R.cos() , cos , 10*eps ) ;        // cos
    EXPECT_NEAR( R.sin() , sin , 10*eps ) ;        // sin
    EXPECT_NEAR( R.theta() , 2*theta , 10*eps ) ;  // theta
  }

  {
    const T theta1 = pi/2 ;
    const qclab::QAngle< T >  angle1( theta1 ) ;
    qclab::QRotation< T >  R1( angle1 ) ;

    const T theta2 = pi/3 ;
    const qclab::QAngle< T >  angle2( theta2 ) ;
    qclab::QRotation< T >  R2( angle2 ) ;

    // operator /
    qclab::QRotation< T > R = R1 / R2 ;
    qclab::QAngle< T >  angle = angle1 - angle2 ;
    T cos   = angle.cos() ;
    T sin   = angle.sin() ;
    T theta = angle.theta() ;
    EXPECT_NEAR( R.cos() , cos , 10*eps ) ;        // cos
    EXPECT_NEAR( R.sin() , sin , 10*eps ) ;        // sin
    EXPECT_NEAR( R.theta() , 2*theta , 10*eps ) ;  // theta
  }

  {
    const T theta1 = pi/2 ;
    const qclab::QAngle< T >  angle1( theta1 ) ;
    qclab::QRotation< T >  R1( angle1 ) ;
    T cos   = R1.cos() ;
    T sin   = R1.sin() ;
    T theta = R1.theta() ;

    // inv
    qclab::QRotation< T >  R2 = R1.inv() ;
    EXPECT_NEAR( R1.cos() ,  cos , eps ) ;      // cos
    EXPECT_NEAR( R2.cos() ,  cos , eps ) ;      // cos
    EXPECT_NEAR( R1.sin() ,  sin , eps ) ;      // sin
    EXPECT_NEAR( R2.sin() , -sin , eps ) ;      // sin
    EXPECT_NEAR( R1.theta() ,  theta , eps ) ;  // theta
    EXPECT_NEAR( R2.theta() , -theta , eps ) ;  // theta
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

