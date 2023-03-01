#include <gtest/gtest.h>
#include "qclab/QAngle.hpp"

template <typename T>
void test_qclab_QAngle() {

  const T pi = 4 * std::atan(1) ;
  const T tol = 10 * std::numeric_limits< T >::epsilon() ;

  {
    qclab::QAngle< T >  angle ;

    EXPECT_EQ( angle.cos() , 1 ) ;    // cos
    EXPECT_EQ( angle.sin() , 0 ) ;    // sin
    EXPECT_EQ( angle.theta() , 0 ) ;  // theta

    // update(angle)
    qclab::QAngle< T >  new_angle( 1 ) ;
    angle.update( new_angle ) ;
    EXPECT_NEAR( angle.cos() , std::cos( 1 ) , tol ) ;
    EXPECT_NEAR( angle.sin() , std::sin( 1 ) , tol ) ;
    EXPECT_NEAR( angle.theta() , 1 , tol ) ;

    // update(theta)
    angle.update( pi/2 ) ;
    EXPECT_EQ( angle.cos() , std::cos( pi/2 ) ) ;
    EXPECT_EQ( angle.sin() , std::sin( pi/2 ) ) ;
    EXPECT_NEAR( angle.theta() , pi/2 , tol ) ;

    // update(cos,sin)
    angle.update( std::cos( pi/3 ) , std::sin( pi/3 ) ) ;
    EXPECT_EQ( angle.cos() , std::cos( pi/3 ) ) ;
    EXPECT_EQ( angle.sin() , std::sin( pi/3 ) ) ;
    EXPECT_NEAR( angle.theta() , pi/3 , tol ) ;

    // operators == and !=
    EXPECT_TRUE( angle != new_angle ) ;
    EXPECT_FALSE( angle == new_angle ) ;
    new_angle.update( angle ) ;
    EXPECT_TRUE( angle == new_angle ) ;
    EXPECT_FALSE( angle != new_angle ) ;
  }

  {
    const T theta = 1 ;
    const T cos = std::cos( theta ) ;
    const T sin = std::sin( theta ) ;
    qclab::QAngle< T >  angle( theta ) ;

    EXPECT_NEAR( angle.cos() , cos , tol ) ;      // cos
    EXPECT_NEAR( angle.sin() , sin , tol ) ;      // sin
    EXPECT_NEAR( angle.theta() , theta , tol ) ;  // theta
  }

  {
    const T theta = pi/2 ;
    const T cos = std::cos( theta ) ;
    const T sin = std::sin( theta ) ;
    qclab::QAngle< T >  angle( theta ) ;

    EXPECT_NEAR( angle.cos() , cos , tol ) ;      // cos
    EXPECT_NEAR( angle.sin() , sin , tol ) ;      // sin
    EXPECT_NEAR( angle.theta() , theta , tol ) ;  // theta
  }

  {
    const T theta = 0.5 ;
    const T cos = std::cos( theta ) ;
    const T sin = std::sin( theta ) ;
    qclab::QAngle< T >  angle( cos , sin ) ;

    EXPECT_NEAR( angle.cos() , cos , tol ) ;      // cos
    EXPECT_NEAR( angle.sin() , sin , tol ) ;      // sin
    EXPECT_NEAR( angle.theta() , theta , tol ) ;  // theta
  }

  {
    const T theta = pi/4 ;
    const T cos = std::cos( theta ) ;
    const T sin = std::sin( theta ) ;
    qclab::QAngle< T >  angle( cos , sin ) ;

    EXPECT_EQ( angle.cos() , cos ) ;              // cos
    EXPECT_EQ( angle.sin() , sin ) ;              // sin
    EXPECT_NEAR( angle.theta() , theta , tol ) ;  // theta
  }

  {
    const T theta1 = pi/2 ;
    const T cos1 = std::cos( theta1 ) ;
    const T sin1 = std::sin( theta1 ) ;
    qclab::QAngle< T >  angle1( theta1 ) ;

    const T theta2 = pi/3 ;
    const T cos2 = std::cos( theta2 ) ;
    const T sin2 = std::sin( theta2 ) ;
    qclab::QAngle< T >  angle2( theta2 ) ;

    // operator +=
    angle1 += angle2 ;
    T theta = theta1 + theta2 ;
    T cos = std::cos( theta ) ;
    T sin = std::sin( theta ) ;
    EXPECT_NEAR( angle1.cos() , cos , tol ) ;       // cos
    EXPECT_NEAR( angle1.sin() , sin , tol ) ;       // sin
    EXPECT_NEAR( angle1.theta() , theta , tol ) ;   // theta
    EXPECT_NEAR( angle2.cos() , cos2 , tol ) ;      // cos
    EXPECT_NEAR( angle2.sin() , sin2 , tol ) ;      // sin
    EXPECT_NEAR( angle2.theta() , theta2 , tol ) ;  // theta
  }

  {
    const T theta1 = pi/2 ;
    const T cos1 = std::cos( theta1 ) ;
    const T sin1 = std::sin( theta1 ) ;
    qclab::QAngle< T >  angle1( theta1 ) ;

    const T theta2 = pi/3 ;
    const T cos2 = std::cos( theta2 ) ;
    const T sin2 = std::sin( theta2 ) ;
    qclab::QAngle< T >  angle2( theta2 ) ;

    // operator -=
    angle1 -= angle2 ;
    T theta = theta1 - theta2 ;
    T cos = std::cos( theta ) ;
    T sin = std::sin( theta ) ;
    EXPECT_NEAR( angle1.cos() , cos , tol ) ;       // cos
    EXPECT_NEAR( angle1.sin() , sin , tol ) ;       // sin
    EXPECT_NEAR( angle1.theta() , theta , tol ) ;   // theta
    EXPECT_NEAR( angle2.cos() , cos2 , tol ) ;      // cos
    EXPECT_NEAR( angle2.sin() , sin2 , tol ) ;      // sin
    EXPECT_NEAR( angle2.theta() , theta2 , tol ) ;  // theta
  }

  {
    const T theta1 = pi/2 ;
    const T cos1 = std::cos( theta1 ) ;
    const T sin1 = std::sin( theta1 ) ;
    qclab::QAngle< T >  angle1( theta1 ) ;

    const T theta2 = pi/3 ;
    const T cos2 = std::cos( theta2 ) ;
    const T sin2 = std::sin( theta2 ) ;
    qclab::QAngle< T >  angle2( theta2 ) ;

    // operator +
    qclab::QAngle< T > sum = angle1 + angle2 ;
    T theta = theta1 + theta2 ;
    T cos = std::cos( theta ) ;
    T sin = std::sin( theta ) ;
    EXPECT_NEAR( sum.cos() , cos , tol ) ;      // cos
    EXPECT_NEAR( sum.sin() , sin , tol ) ;      // sin
    EXPECT_NEAR( sum.theta() , theta , tol ) ;  // theta

    // operator -
    qclab::QAngle< T > minus = angle1 - angle2 ;
    theta = theta1 - theta2 ;
    cos = std::cos( theta ) ;
    sin = std::sin( theta ) ;
    EXPECT_NEAR( minus.cos() , cos , tol ) ;      // cos
    EXPECT_NEAR( minus.sin() , sin , tol ) ;      // sin
    EXPECT_NEAR( minus.theta() , theta , tol ) ;  // theta
  }

  {
    const T theta = 1 ;
    const T cos = std::cos( theta ) ;
    const T sin = std::sin( theta ) ;
    qclab::QAngle< T >  angle1( theta ) ;
    qclab::QAngle< T >  angle2 ;

    // operator -
    angle2 = -angle1 ;
    EXPECT_NEAR( angle1.cos() ,  cos , tol ) ;      // cos
    EXPECT_NEAR( angle1.sin() , sin , tol ) ;       // sin
    EXPECT_NEAR( angle1.theta() , theta , tol ) ;   // theta
    EXPECT_NEAR( angle2.cos() ,  cos , tol ) ;      // cos
    EXPECT_NEAR( angle2.sin() , -sin , tol ) ;      // sin
    EXPECT_NEAR( angle2.theta() , -theta , tol ) ;  // theta
  }

}


/*
 * float
 */
TEST( qclab_QAngle , float ) {
  test_qclab_QAngle< float >() ;
}

/*
 * double
 */
TEST( qclab_QAngle , double ) {
  test_qclab_QAngle< double >() ;
}

