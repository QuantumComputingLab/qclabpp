#include <gtest/gtest.h>
#include "QAngle.hpp"

template <typename T>
void test_qclab_QAngle() {

  const T pi = 4 * std::atan(1) ;
  const T eps = std::numeric_limits< T >::epsilon() ;

  {
    qclab::QAngle< T >  angle ;

    EXPECT_EQ( angle.cos() , 1 ) ;    // cos
    EXPECT_EQ( angle.sin() , 0 ) ;    // sin
    EXPECT_EQ( angle.theta() , 0 ) ;  // theta

    // update(angle)
    qclab::QAngle< T >  new_angle( 1 ) ;
    angle.update( new_angle ) ;
    EXPECT_NEAR( angle.cos() , std::cos( 1 ) , eps ) ;
    EXPECT_NEAR( angle.sin() , std::sin( 1 ) , eps ) ;
    EXPECT_NEAR( angle.theta() , 1 , eps ) ;

    // update(theta)
    angle.update( pi/2 ) ;
    EXPECT_EQ( angle.cos() , std::cos( pi/2 ) ) ;
    EXPECT_EQ( angle.sin() , std::sin( pi/2 ) ) ;
    EXPECT_NEAR( angle.theta() , pi/2 , eps ) ;

    // update(cos,sin)
    angle.update( std::cos( pi/3 ) , std::sin( pi/3 ) ) ;
    EXPECT_EQ( angle.cos() , std::cos( pi/3 ) ) ;
    EXPECT_EQ( angle.sin() , std::sin( pi/3 ) ) ;
    EXPECT_NEAR( angle.theta() , pi/3 , eps ) ;

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

    EXPECT_NEAR( angle.cos() , cos , eps ) ;      // cos
    EXPECT_NEAR( angle.sin() , sin , eps ) ;      // sin
    EXPECT_NEAR( angle.theta() , theta , eps ) ;  // theta
  }

  {
    const T theta = pi/2 ;
    const T cos = std::cos( theta ) ;
    const T sin = std::sin( theta ) ;
    qclab::QAngle< T >  angle( theta ) ;

    EXPECT_NEAR( angle.cos() , cos , eps ) ;      // cos
    EXPECT_NEAR( angle.sin() , sin , eps ) ;      // sin
    EXPECT_NEAR( angle.theta() , theta , eps ) ;  // theta
  }

  {
    const T theta = 0.5 ;
    const T cos = std::cos( theta ) ;
    const T sin = std::sin( theta ) ;
    qclab::QAngle< T >  angle( cos , sin ) ;

    EXPECT_NEAR( angle.cos() , cos , eps ) ;      // cos
    EXPECT_NEAR( angle.sin() , sin , eps ) ;      // sin
    EXPECT_NEAR( angle.theta() , theta , eps ) ;  // theta
  }

  {
    const T theta = pi/4 ;
    const T cos = std::cos( theta ) ;
    const T sin = std::sin( theta ) ;
    qclab::QAngle< T >  angle( cos , sin ) ;

    EXPECT_EQ( angle.cos() , cos ) ;              // cos
    EXPECT_EQ( angle.sin() , sin ) ;              // sin
    EXPECT_NEAR( angle.theta() , theta , eps ) ;  // theta
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

