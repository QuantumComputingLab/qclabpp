#include <gtest/gtest.h>
#include "qclab/qgates/RotationZ.hpp"

template <typename T>
void test_qclab_qgates_RotationZ() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  {
    qclab::qgates::RotationZ< T >  Rz ;

    EXPECT_EQ( Rz.nbQubits() , 1 ) ;   // nbQubits
    EXPECT_FALSE( Rz.fixed() ) ;       // fixed
    EXPECT_FALSE( Rz.controlled() ) ;  // controlled
    EXPECT_EQ( Rz.cos() , 1.0 ) ;      // cos
    EXPECT_EQ( Rz.sin() , 0.0 ) ;      // sin
    EXPECT_EQ( Rz.theta() , 0.0 ) ;    // theta

    // matrix
    auto eye = qclab::dense::eye< T >( 2 ) ;
    EXPECT_TRUE( Rz.matrix() == eye ) ;

    // qubit
    EXPECT_EQ( Rz.qubit() , 0 ) ;
    Rz.setQubit( 2 ) ;
    EXPECT_EQ( Rz.qubit() , 2 ) ;

    // qubits
    auto qubits = Rz.qubits() ;
    EXPECT_EQ( qubits.size() , 1 ) ;
    EXPECT_EQ( qubits[0] , 2 ) ;
    int qnew = 3 ;
    Rz.setQubits( &qnew ) ;
    EXPECT_EQ( Rz.qubit() , 3 ) ;

    // fixed
    Rz.makeFixed() ;
    EXPECT_TRUE( Rz.fixed() ) ;
    Rz.makeVariable() ;
    EXPECT_FALSE( Rz.fixed() ) ;

    // update(rot)
    qclab::QRotation< R >  new_rot( 1 ) ;
    Rz.update( new_rot ) ;
    EXPECT_NEAR( Rz.theta() , 1 , eps ) ;
    EXPECT_NEAR( Rz.cos() , std::cos( 0.5 ) , eps ) ;
    EXPECT_NEAR( Rz.sin() , std::sin( 0.5 ) , eps ) ;

    // update(theta)
    Rz.update( pi/2 ) ;
    EXPECT_NEAR( Rz.theta() , pi/2 , eps ) ;
    EXPECT_NEAR( Rz.cos() , std::cos( pi/4 ) , eps ) ;
    EXPECT_NEAR( Rz.sin() , std::sin( pi/4 ) , eps ) ;

    // matrix
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    EXPECT_NEAR( std::real( Rz.matrix()(0,0) ) ,  cos , eps ) ;
    EXPECT_NEAR( std::imag( Rz.matrix()(0,0) ) , -sin , eps ) ;
    EXPECT_EQ( Rz.matrix()(1,0) , T(0.0) ) ;
    EXPECT_EQ( Rz.matrix()(0,1) , T(0.0) ) ;
    EXPECT_NEAR( std::real( Rz.matrix()(1,1) ) ,  cos , eps ) ;
    EXPECT_NEAR( std::imag( Rz.matrix()(1,1) ) ,  sin , eps ) ;

    // print
    Rz.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( Rz.toQASM( qasm ) , 0 ) ;
    std::string qasm_check = "rz(" + qclab::qasm( Rz.theta() ) + ") q[3];\n" ;
    EXPECT_EQ( qasm.str() , qasm_check ) ;
    std::cout << qasm.str() ;

    // update(cos,sin)
    Rz.update( std::cos( pi/3 ) , std::sin( pi/3 ) ) ;
    EXPECT_EQ( Rz.cos() , std::cos( pi/3 ) ) ;
    EXPECT_EQ( Rz.sin() , std::sin( pi/3 ) ) ;
    EXPECT_NEAR( Rz.theta() , 2*(pi/3) , eps ) ;

    // operators == and !=
    qclab::qgates::RotationZ< T >  Rz2( std::cos( pi/3 ) , std::sin( pi/3 ) ) ;
    EXPECT_TRUE( Rz == Rz2 ) ;
    EXPECT_FALSE( Rz != Rz2 ) ;
    Rz2.update( 1 ) ;
    EXPECT_TRUE( Rz != Rz2 ) ;
    EXPECT_FALSE( Rz == Rz2 ) ;
  }

  {
    qclab::qgates::RotationZ< T >  Rz( pi/2 ) ;

    EXPECT_EQ( Rz.nbQubits() , 1 ) ;                    // nbQubits
    EXPECT_FALSE( Rz.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rz.controlled() ) ;                   // controlled
    EXPECT_EQ( Rz.qubit() , 0 ) ;                       // qubit

    EXPECT_NEAR( Rz.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Rz.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( Rz.sin() , std::sin( pi/4 ) , eps ) ;  // sin
  }

  {
    qclab::qgates::RotationZ< T >  Rz( 5 , pi/2 ) ;

    EXPECT_EQ( Rz.nbQubits() , 1 ) ;                    // nbQubits
    EXPECT_FALSE( Rz.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rz.controlled() ) ;                   // controlled
    EXPECT_EQ( Rz.qubit() , 5 ) ;                       // qubit

    EXPECT_NEAR( Rz.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Rz.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( Rz.sin() , std::sin( pi/4 ) , eps ) ;  // sin
  }

  {
    qclab::qgates::RotationZ< T >  Rz( 5 , pi/2 , true ) ;

    EXPECT_EQ( Rz.nbQubits() , 1 ) ;                    // nbQubits
    EXPECT_TRUE( Rz.fixed() ) ;                         // fixed
    EXPECT_FALSE( Rz.controlled() ) ;                   // controlled
    EXPECT_EQ( Rz.qubit() , 5 ) ;                       // qubit

    EXPECT_NEAR( Rz.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Rz.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( Rz.sin() , std::sin( pi/4 ) , eps ) ;  // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationZ< T >  Rz( cos , sin ) ;

    EXPECT_EQ( Rz.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_FALSE( Rz.fixed() ) ;              // fixed
    EXPECT_FALSE( Rz.controlled() ) ;         // controlled
    EXPECT_EQ( Rz.qubit() , 0 ) ;             // qubit

    EXPECT_NEAR( Rz.theta() , pi/2 , eps ) ;  // theta
    EXPECT_NEAR( Rz.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( Rz.sin() , sin , eps ) ;     // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationZ< T >  Rz( 5 , cos , sin ) ;

    EXPECT_EQ( Rz.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_FALSE( Rz.fixed() ) ;              // fixed
    EXPECT_FALSE( Rz.controlled() ) ;         // controlled
    EXPECT_EQ( Rz.qubit() , 5 ) ;             // qubit

    EXPECT_NEAR( Rz.theta() , pi/2 , eps ) ;  // theta
    EXPECT_NEAR( Rz.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( Rz.sin() , sin , eps ) ;     // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationZ< T >  Rz( 5 , cos , sin , true ) ;

    EXPECT_EQ( Rz.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_TRUE( Rz.fixed() ) ;               // fixed
    EXPECT_FALSE( Rz.controlled() ) ;         // controlled
    EXPECT_EQ( Rz.qubit() , 5 ) ;             // qubit

    EXPECT_NEAR( Rz.theta() , pi/2 , eps ) ;  // theta
    EXPECT_NEAR( Rz.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( Rz.sin() , sin , eps ) ;     // sin
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >       angle1( theta1 ) ;
    const qclab::QRotation< R >    rot1( angle1 ) ;
    qclab::qgates::RotationZ< T >  R1( rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >       angle2( theta2 ) ;
    const qclab::QRotation< R >    rot2( angle2 ) ;
    qclab::qgates::RotationZ< T >  R2( rot2 ) ;

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
    const qclab::QAngle< R >       angle1( theta1 ) ;
    const qclab::QRotation< R >    rot1( angle1 ) ;
    qclab::qgates::RotationZ< T >  R1( rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >       angle2( theta2 ) ;
    const qclab::QRotation< R >    rot2( angle2 ) ;
    qclab::qgates::RotationZ< T >  R2( rot2 ) ;

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
    const qclab::QAngle< R >       angle1( theta1 ) ;
    const qclab::QRotation< R >    rot1( angle1 ) ;
    qclab::qgates::RotationZ< T >  R1( rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >       angle2( theta2 ) ;
    const qclab::QRotation< R >    rot2( angle2 ) ;
    qclab::qgates::RotationZ< T >  R2( rot2 ) ;

    // operator *
    qclab::qgates::RotationZ< T > R12 = R1 * R2 ;
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
    const qclab::QAngle< R >       angle1( theta1 ) ;
    const qclab::QRotation< R >    rot1( angle1 ) ;
    qclab::qgates::RotationZ< T >  R1( rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >       angle2( theta2 ) ;
    const qclab::QRotation< R >    rot2( angle2 ) ;
    qclab::qgates::RotationZ< T >  R2( rot2 ) ;

    // operator /
    qclab::qgates::RotationZ< T > R12 = R1 / R2 ;
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
    const qclab::QAngle< R >       angle1( theta1 ) ;
    const qclab::QRotation< R >    rot1( angle1 ) ;
    qclab::qgates::RotationZ< T >  R1( rot1 ) ;
    R cos   = R1.cos() ;
    R sin   = R1.sin() ;
    R theta = R1.theta() ;

    // inv
    qclab::qgates::RotationZ< T >  R2 = R1.inv() ;
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
TEST( qclab_qgates_RotationZ , complex_float ) {
  test_qclab_qgates_RotationZ< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_RotationZ , complex_double ) {
  test_qclab_qgates_RotationZ< std::complex< double > >() ;
}

