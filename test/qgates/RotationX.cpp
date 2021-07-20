#include <gtest/gtest.h>
#include "qclab/qgates/RotationX.hpp"

template <typename T>
void test_qclab_qgates_RotationX() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  {
    qclab::qgates::RotationX< T >  Rx ;

    EXPECT_EQ( Rx.nbQubits() , 1 ) ;   // nbQubits
    EXPECT_FALSE( Rx.fixed() ) ;       // fixed
    EXPECT_FALSE( Rx.controlled() ) ;  // controlled
    EXPECT_EQ( Rx.cos() , 1.0 ) ;      // cos
    EXPECT_EQ( Rx.sin() , 0.0 ) ;      // sin
    EXPECT_EQ( Rx.theta() , 0.0 ) ;    // theta

    // matrix
    auto eye = qclab::dense::eye< T >( 2 ) ;
    EXPECT_TRUE( Rx.matrix() == eye ) ;

    // qubit
    EXPECT_EQ( Rx.qubit() , 0 ) ;
    Rx.setQubit( 2 ) ;
    EXPECT_EQ( Rx.qubit() , 2 ) ;

    // qubits
    auto qubits = Rx.qubits() ;
    EXPECT_EQ( qubits.size() , 1 ) ;
    EXPECT_EQ( qubits[0] , 2 ) ;
    int qnew = 3 ;
    Rx.setQubits( &qnew ) ;
    EXPECT_EQ( Rx.qubit() , 3 ) ;

    // fixed
    Rx.makeFixed() ;
    EXPECT_TRUE( Rx.fixed() ) ;
    Rx.makeVariable() ;
    EXPECT_FALSE( Rx.fixed() ) ;

    // update(angle)
    qclab::QAngle< R >  new_angle( 0.5 ) ;
    Rx.update( new_angle ) ;
    EXPECT_NEAR( Rx.theta() , 1 , eps ) ;
    EXPECT_NEAR( Rx.cos() , std::cos( 0.5 ) , eps ) ;
    EXPECT_NEAR( Rx.sin() , std::sin( 0.5 ) , eps ) ;

    // update(theta)
    Rx.update( pi/2 ) ;
    EXPECT_NEAR( Rx.theta() , pi/2 , eps ) ;
    EXPECT_NEAR( Rx.cos() , std::cos( pi/4 ) , eps ) ;
    EXPECT_NEAR( Rx.sin() , std::sin( pi/4 ) , eps ) ;

    // matrix
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    EXPECT_EQ( Rx.matrix()(0,0) , T(cos,0.0) ) ;
    EXPECT_EQ( Rx.matrix()(1,0) , T(0.0,-sin) ) ;
    EXPECT_EQ( Rx.matrix()(0,1) , T(0.0,-sin) ) ;
    EXPECT_EQ( Rx.matrix()(1,1) , T(cos,0.0) ) ;

    // print
    Rx.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( Rx.toQASM( qasm ) , 0 ) ;
    std::string qasm_check = "rx(" + qclab::qasm( Rx.theta() ) + ") q[3];\n" ;
    EXPECT_EQ( qasm.str() , qasm_check ) ;
    std::cout << qasm.str() ;

    // update(cos,sin)
    Rx.update( std::cos( pi/3 ) , std::sin( pi/3 ) ) ;
    EXPECT_NEAR( Rx.cos() , std::cos( pi/3 ) , eps ) ;
    EXPECT_NEAR( Rx.sin() , std::sin( pi/3 ) , eps ) ;
    EXPECT_NEAR( Rx.theta() , 2*(pi/3) , eps ) ;

    // operators == and !=
    qclab::qgates::RotationX< T >  Rx2( std::cos( pi/3 ) , std::sin( pi/3 ) ) ;
    EXPECT_TRUE( Rx == Rx2 ) ;
    EXPECT_FALSE( Rx != Rx2 ) ;
    Rx2.update( 1 ) ;
    EXPECT_TRUE( Rx != Rx2 ) ;
    EXPECT_FALSE( Rx == Rx2 ) ;
  }

  {
    qclab::qgates::RotationX< T >  Rx( pi/2 ) ;

    EXPECT_EQ( Rx.nbQubits() , 1 ) ;                    // nbQubits
    EXPECT_FALSE( Rx.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rx.controlled() ) ;                   // controlled
    EXPECT_EQ( Rx.qubit() , 0 ) ;                       // qubit

    EXPECT_NEAR( Rx.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Rx.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( Rx.sin() , std::sin( pi/4 ) , eps ) ;  // sin
  }

  {
    qclab::qgates::RotationX< T >  Rx( 5 , pi/2 ) ;

    EXPECT_EQ( Rx.nbQubits() , 1 ) ;                    // nbQubits
    EXPECT_FALSE( Rx.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rx.controlled() ) ;                   // controlled
    EXPECT_EQ( Rx.qubit() , 5 ) ;                       // qubit

    EXPECT_NEAR( Rx.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Rx.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( Rx.sin() , std::sin( pi/4 ) , eps ) ;  // sin
  }

  {
    qclab::qgates::RotationX< T >  Rx( 5 , pi/2 , true ) ;

    EXPECT_EQ( Rx.nbQubits() , 1 ) ;                    // nbQubits
    EXPECT_TRUE( Rx.fixed() ) ;                         // fixed
    EXPECT_FALSE( Rx.controlled() ) ;                   // controlled
    EXPECT_EQ( Rx.qubit() , 5 ) ;                       // qubit

    EXPECT_NEAR( Rx.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Rx.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( Rx.sin() , std::sin( pi/4 ) , eps ) ;  // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationX< T >  Rx( cos , sin ) ;

    EXPECT_EQ( Rx.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_FALSE( Rx.fixed() ) ;              // fixed
    EXPECT_FALSE( Rx.controlled() ) ;         // controlled
    EXPECT_EQ( Rx.qubit() , 0 ) ;             // qubit

    EXPECT_NEAR( Rx.theta() , pi/2 , eps ) ;  // theta
    EXPECT_NEAR( Rx.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( Rx.sin() , sin , eps ) ;     // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationX< T >  Rx( 5 , cos , sin ) ;

    EXPECT_EQ( Rx.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_FALSE( Rx.fixed() ) ;              // fixed
    EXPECT_FALSE( Rx.controlled() ) ;         // controlled
    EXPECT_EQ( Rx.qubit() , 5 ) ;             // qubit

    EXPECT_NEAR( Rx.theta() , pi/2 , eps ) ;  // theta
    EXPECT_NEAR( Rx.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( Rx.sin() , sin , eps ) ;     // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationX< T >  Rx( 5 , cos , sin , true ) ;

    EXPECT_EQ( Rx.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_TRUE( Rx.fixed() ) ;               // fixed
    EXPECT_FALSE( Rx.controlled() ) ;         // controlled
    EXPECT_EQ( Rx.qubit() , 5 ) ;             // qubit

    EXPECT_NEAR( Rx.theta() , pi/2 , eps ) ;  // theta
    EXPECT_NEAR( Rx.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( Rx.sin() , sin , eps ) ;     // sin
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >  angle1( theta1 ) ;
    qclab::qgates::RotationX< T >  R1( angle1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >  angle2( theta2 ) ;
    qclab::qgates::RotationX< T >  R2( angle2 ) ;

    // operator *=
    R1 *= R2 ;
    qclab::QAngle< R >  angle = angle1 + angle2 ;
    R cos   = angle.cos() ;
    R sin   = angle.sin() ;
    R theta = angle.theta() ;
    EXPECT_NEAR( R1.cos() , cos , 10*eps ) ;                 // cos
    EXPECT_NEAR( R1.sin() , sin , 10*eps ) ;                 // sin
    EXPECT_NEAR( R1.theta() , 2*theta , 10*eps ) ;           // theta
    EXPECT_NEAR( R2.cos() , std::cos( theta2 ) , 10*eps ) ;  // cos
    EXPECT_NEAR( R2.sin() , std::sin( theta2 ) , 10*eps ) ;  // sin
    EXPECT_NEAR( R2.theta() , 2*theta2 , 10*eps ) ;          // theta
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >  angle1( theta1 ) ;
    qclab::qgates::RotationX< T >  R1( angle1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >  angle2( theta2 ) ;
    qclab::qgates::RotationX< T >  R2( angle2 ) ;

    // operator /=
    R1 /= R2 ;
    qclab::QAngle< R >  angle = angle1 - angle2 ;
    R cos   = angle.cos() ;
    R sin   = angle.sin() ;
    R theta = angle.theta() ;
    EXPECT_NEAR( R1.cos() , cos , 10*eps ) ;                 // cos
    EXPECT_NEAR( R1.sin() , sin , 10*eps ) ;                 // sin
    EXPECT_NEAR( R1.theta() , 2*theta , 10*eps ) ;           // theta
    EXPECT_NEAR( R2.cos() , std::cos( theta2 ) , 10*eps ) ;  // cos
    EXPECT_NEAR( R2.sin() , std::sin( theta2 ) , 10*eps ) ;  // sin
    EXPECT_NEAR( R2.theta() , 2*theta2 , 10*eps ) ;          // theta
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >  angle1( theta1 ) ;
    qclab::qgates::RotationX< T >  R1( angle1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >  angle2( theta2 ) ;
    qclab::qgates::RotationX< T >  R2( angle2 ) ;

    // operator *
    qclab::qgates::RotationX< T > R12 = R1 * R2 ;
    qclab::QAngle< R >  angle = angle1 + angle2 ;
    R cos   = angle.cos() ;
    R sin   = angle.sin() ;
    R theta = angle.theta() ;
    EXPECT_NEAR( R12.cos() , cos , 10*eps ) ;        // cos
    EXPECT_NEAR( R12.sin() , sin , 10*eps ) ;        // sin
    EXPECT_NEAR( R12.theta() , 2*theta , 10*eps ) ;  // theta
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >  angle1( theta1 ) ;
    qclab::qgates::RotationX< T >  R1( angle1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >  angle2( theta2 ) ;
    qclab::qgates::RotationX< T >  R2( angle2 ) ;

    // operator /
    qclab::qgates::RotationX< T > R12 = R1 / R2 ;
    qclab::QAngle< R >  angle = angle1 - angle2 ;
    R cos   = angle.cos() ;
    R sin   = angle.sin() ;
    R theta = angle.theta() ;
    EXPECT_NEAR( R12.cos() , cos , 10*eps ) ;        // cos
    EXPECT_NEAR( R12.sin() , sin , 10*eps ) ;        // sin
    EXPECT_NEAR( R12.theta() , 2*theta , 10*eps ) ;  // theta
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >  angle1( theta1 ) ;
    qclab::qgates::RotationX< T >  R1( angle1 ) ;
    R cos   = R1.cos() ;
    R sin   = R1.sin() ;
    R theta = R1.theta() ;

    // inv
    qclab::qgates::RotationX< T >  R2 = R1.inv() ;
    EXPECT_NEAR( R1.cos() ,  cos , eps ) ;      // cos
    EXPECT_NEAR( R2.cos() ,  cos , eps ) ;      // cos
    EXPECT_NEAR( R1.sin() ,  sin , eps ) ;      // sin
    EXPECT_NEAR( R2.sin() , -sin , eps ) ;      // sin
    EXPECT_NEAR( R1.theta() ,  theta , eps ) ;  // theta
    EXPECT_NEAR( R2.theta() , -theta , eps ) ;  // theta
  }

}


/*
 * complex float
 */
TEST( qclab_qgates_RotationX , complex_float ) {
  test_qclab_qgates_RotationX< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_RotationX , complex_double ) {
  test_qclab_qgates_RotationX< std::complex< double > >() ;
}

