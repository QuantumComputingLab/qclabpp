#include <gtest/gtest.h>
#include "qclab/qgates/RotationZZ.hpp"

template <typename T>
void test_qclab_qgates_RotationZZ() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R tol = 10 * std::numeric_limits< R >::epsilon() ;

  {
    qclab::qgates::RotationZZ< T >  Rzz ;

    EXPECT_EQ( Rzz.nbQubits() , 2 ) ;   // nbQubits
    EXPECT_FALSE( Rzz.fixed() ) ;       // fixed
    EXPECT_FALSE( Rzz.controlled() ) ;  // controlled
    EXPECT_EQ( Rzz.cos() , 1.0 ) ;      // cos
    EXPECT_EQ( Rzz.sin() , 0.0 ) ;      // sin
    EXPECT_EQ( Rzz.theta() , 0.0 ) ;    // theta

    // matrix
    auto eye = qclab::dense::eye< T >( 4 ) ;
    EXPECT_TRUE( Rzz.matrix() == eye ) ;

    // qubits
    EXPECT_EQ( Rzz.qubits().size() , 2 ) ;
    EXPECT_EQ( Rzz.qubits()[0] , 0 ) ;
    EXPECT_EQ( Rzz.qubits()[1] , 1 ) ;
    int qnew[2] = { 3 , 5 } ;
    Rzz.setQubits( &qnew[0] ) ;
    EXPECT_EQ( Rzz.qubits().size() , 2 ) ;
    EXPECT_EQ( Rzz.qubits()[0] , 3 ) ;
    EXPECT_EQ( Rzz.qubits()[1] , 5 ) ;

    // fixed
    Rzz.makeFixed() ;
    EXPECT_TRUE( Rzz.fixed() ) ;
    Rzz.makeVariable() ;
    EXPECT_FALSE( Rzz.fixed() ) ;

    // update(rot)
    qclab::QRotation< R >  new_rot( 1 ) ;
    Rzz.update( new_rot ) ;
    EXPECT_NEAR( Rzz.theta() , 1 , tol ) ;
    EXPECT_NEAR( Rzz.cos() , std::cos( 0.5 ) , tol ) ;
    EXPECT_NEAR( Rzz.sin() , std::sin( 0.5 ) , tol ) ;

    // update(theta)
    Rzz.update( pi/2 ) ;
    EXPECT_NEAR( Rzz.theta() , pi/2 , tol ) ;
    EXPECT_NEAR( Rzz.cos() , std::cos( pi/4 ) , tol ) ;
    EXPECT_NEAR( Rzz.sin() , std::sin( pi/4 ) , tol ) ;

    // matrix
    const T a = std::exp( T(0,-pi/4) ) ;
    const T b = std::exp( T(0, pi/4) ) ;
    EXPECT_NEAR( std::real( Rzz.matrix()(0,0) ) , std::real( a ) , tol ) ;
    EXPECT_NEAR( std::imag( Rzz.matrix()(0,0) ) , std::imag( a ) , tol ) ;
    EXPECT_NEAR( std::real( Rzz.matrix()(1,1) ) , std::real( b ) , tol ) ;
    EXPECT_NEAR( std::imag( Rzz.matrix()(1,1) ) , std::imag( b ) , tol ) ;
    EXPECT_NEAR( std::real( Rzz.matrix()(2,2) ) , std::real( b ) , tol ) ;
    EXPECT_NEAR( std::imag( Rzz.matrix()(2,2) ) , std::imag( b ) , tol ) ;
    EXPECT_NEAR( std::real( Rzz.matrix()(3,3) ) , std::real( a ) , tol ) ;
    EXPECT_NEAR( std::imag( Rzz.matrix()(3,3) ) , std::imag( a ) , tol ) ;
    EXPECT_EQ( Rzz.matrix()(1,0) , T(0) ) ;
    EXPECT_EQ( Rzz.matrix()(2,0) , T(0) ) ;
    EXPECT_EQ( Rzz.matrix()(3,0) , T(0) ) ;
    EXPECT_EQ( Rzz.matrix()(0,1) , T(0) ) ;
    EXPECT_EQ( Rzz.matrix()(2,1) , T(0) ) ;
    EXPECT_EQ( Rzz.matrix()(3,1) , T(0) ) ;
    EXPECT_EQ( Rzz.matrix()(0,2) , T(0) ) ;
    EXPECT_EQ( Rzz.matrix()(1,2) , T(0) ) ;
    EXPECT_EQ( Rzz.matrix()(3,2) , T(0) ) ;
    EXPECT_EQ( Rzz.matrix()(0,3) , T(0) ) ;
    EXPECT_EQ( Rzz.matrix()(1,3) , T(0) ) ;
    EXPECT_EQ( Rzz.matrix()(2,3) , T(0) ) ;

    // print
    Rzz.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( Rzz.toQASM( qasm ) , 0 ) ;
    auto str_theta = qclab::qasm( Rzz.theta() ) ;
    std::string qasm_check = "rzz(" + str_theta + ") q[3], q[5];\n" ;
    EXPECT_EQ( qasm.str() , qasm_check ) ;
    std::cout << qasm.str() ;

    // update(cos,sin)
    Rzz.update( std::cos( pi/3 ) , std::sin( pi/3 ) ) ;
    EXPECT_NEAR( Rzz.cos() , std::cos( pi/3 ) , tol ) ;
    EXPECT_NEAR( Rzz.sin() , std::sin( pi/3 ) , tol ) ;
    EXPECT_NEAR( Rzz.theta() , 2*(pi/3) , tol ) ;

    // operators == and !=
    qclab::qgates::RotationZZ< T >  Rzz2( std::cos( pi/3 ) , std::sin( pi/3 ) );
    EXPECT_TRUE( Rzz == Rzz2 ) ;
    EXPECT_FALSE( Rzz != Rzz2 ) ;
    Rzz2.update( 1 ) ;
    EXPECT_TRUE( Rzz != Rzz2 ) ;
    EXPECT_FALSE( Rzz == Rzz2 ) ;
  }

  //
  // constructors without qubits
  //
  {
    qclab::QRotation< R >           rot( pi/2 ) ;
    qclab::qgates::RotationZZ< T >  Rzz( rot ) ;

    EXPECT_EQ( Rzz.nbQubits() , 2 ) ;                    // nbQubits
    EXPECT_FALSE( Rzz.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rzz.controlled() ) ;                   // controlled

    auto qubits = Rzz.qubits() ;
    EXPECT_EQ( qubits[0] , 0 ) ;                         // qubit0
    EXPECT_EQ( qubits[1] , 1 ) ;                         // qubit1

    EXPECT_NEAR( Rzz.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Rzz.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Rzz.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    qclab::qgates::RotationZZ< T >  Rzz( pi/2 ) ;

    EXPECT_EQ( Rzz.nbQubits() , 2 ) ;                    // nbQubits
    EXPECT_FALSE( Rzz.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rzz.controlled() ) ;                   // controlled

    auto qubits = Rzz.qubits() ;
    EXPECT_EQ( qubits[0] , 0 ) ;                         // qubit0
    EXPECT_EQ( qubits[1] , 1 ) ;                         // qubit1

    EXPECT_NEAR( Rzz.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Rzz.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Rzz.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationZZ< T >  Rzz( cos , sin ) ;

    EXPECT_EQ( Rzz.nbQubits() , 2 ) ;          // nbQubits
    EXPECT_FALSE( Rzz.fixed() ) ;              // fixed
    EXPECT_FALSE( Rzz.controlled() ) ;         // controlled

    auto qubits = Rzz.qubits() ;
    EXPECT_EQ( qubits[0] , 0 ) ;               // qubit0
    EXPECT_EQ( qubits[1] , 1 ) ;               // qubit1

    EXPECT_NEAR( Rzz.theta() , pi/2 , tol ) ;  // theta
    EXPECT_NEAR( Rzz.cos() , cos , tol ) ;     // cos
    EXPECT_NEAR( Rzz.sin() , sin , tol ) ;     // sin
  }


  //
  // constructors with *qubits
  //
  {
    int qbits[2] = { 3 , 5 } ;
    qclab::QRotation< R >           rot( pi/2 ) ;
    qclab::qgates::RotationZZ< T >  Rzz( &qbits[0] , rot ) ;

    EXPECT_EQ( Rzz.nbQubits() , 2 ) ;                    // nbQubits
    EXPECT_FALSE( Rzz.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rzz.controlled() ) ;                   // controlled

    auto qubits = Rzz.qubits() ;
    EXPECT_EQ( qubits[0] , 3 ) ;                         // qubit0
    EXPECT_EQ( qubits[1] , 5 ) ;                         // qubit1

    EXPECT_NEAR( Rzz.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Rzz.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Rzz.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    int qbits[2] = { 3 , 5 } ;
    qclab::qgates::RotationZZ< T >  Rzz( &qbits[0] , pi/2 ) ;

    EXPECT_EQ( Rzz.nbQubits() , 2 ) ;                    // nbQubits
    EXPECT_FALSE( Rzz.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rzz.controlled() ) ;                   // controlled

    auto qubits = Rzz.qubits() ;
    EXPECT_EQ( qubits[0] , 3 ) ;                         // qubit0
    EXPECT_EQ( qubits[1] , 5 ) ;                         // qubit1

    EXPECT_NEAR( Rzz.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Rzz.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Rzz.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    int qbits[2] = { 3 , 5 } ;
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationZZ< T >  Rzz( &qbits[0] , cos , sin ) ;

    EXPECT_EQ( Rzz.nbQubits() , 2 ) ;          // nbQubits
    EXPECT_FALSE( Rzz.fixed() ) ;              // fixed
    EXPECT_FALSE( Rzz.controlled() ) ;         // controlled

    auto qubits = Rzz.qubits() ;
    EXPECT_EQ( qubits[0] , 3 ) ;               // qubit0
    EXPECT_EQ( qubits[1] , 5 ) ;               // qubit1

    EXPECT_NEAR( Rzz.theta() , pi/2 , tol ) ;  // theta
    EXPECT_NEAR( Rzz.cos() , cos , tol ) ;     // cos
    EXPECT_NEAR( Rzz.sin() , sin , tol ) ;     // sin
  }


  //
  // constructors with qubit0 and qubit1
  //
  {
    int qubit0 = 3 ;
    int qubit1 = 5 ;
    qclab::QRotation< R >           rot( pi/2 ) ;
    qclab::qgates::RotationZZ< T >  Rzz( qubit0 , qubit1 , rot ) ;

    EXPECT_EQ( Rzz.nbQubits() , 2 ) ;                    // nbQubits
    EXPECT_FALSE( Rzz.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rzz.controlled() ) ;                   // controlled

    auto qubits = Rzz.qubits() ;
    EXPECT_EQ( qubits[0] , 3 ) ;                         // qubit0
    EXPECT_EQ( qubits[1] , 5 ) ;                         // qubit1

    EXPECT_NEAR( Rzz.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Rzz.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Rzz.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    int qubit0 = 3 ;
    int qubit1 = 5 ;
    qclab::qgates::RotationZZ< T >  Rzz( qubit0 , qubit1 , pi/2 ) ;

    EXPECT_EQ( Rzz.nbQubits() , 2 ) ;                    // nbQubits
    EXPECT_FALSE( Rzz.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rzz.controlled() ) ;                   // controlled

    auto qubits = Rzz.qubits() ;
    EXPECT_EQ( qubits[0] , 3 ) ;                         // qubit0
    EXPECT_EQ( qubits[1] , 5 ) ;                         // qubit1

    EXPECT_NEAR( Rzz.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Rzz.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Rzz.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    int qubit0 = 3 ;
    int qubit1 = 5 ;
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationZZ< T >  Rzz( qubit0 , qubit1 , cos , sin ) ;

    EXPECT_EQ( Rzz.nbQubits() , 2 ) ;          // nbQubits
    EXPECT_FALSE( Rzz.fixed() ) ;              // fixed
    EXPECT_FALSE( Rzz.controlled() ) ;         // controlled

    auto qubits = Rzz.qubits() ;
    EXPECT_EQ( qubits[0] , 3 ) ;               // qubit0
    EXPECT_EQ( qubits[1] , 5 ) ;               // qubit1

    EXPECT_NEAR( Rzz.theta() , pi/2 , tol ) ;  // theta
    EXPECT_NEAR( Rzz.cos() , cos , tol ) ;     // cos
    EXPECT_NEAR( Rzz.sin() , sin , tol ) ;     // sin
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >        angle1( theta1 ) ;
    const qclab::QRotation< R >     rot1( angle1 ) ;
    qclab::qgates::RotationZZ< T >  R1( 0 , 1 , rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >        angle2( theta2 ) ;
    const qclab::QRotation< R >     rot2( angle2 ) ;
    qclab::qgates::RotationZZ< T >  R2( 0 , 1 , rot2 ) ;

    // operator *=
    R1 *= R2 ;
    qclab::QAngle< R >  angle = angle1 + angle2 ;
    R cos   = angle.cos() ;
    R sin   = angle.sin() ;
    R theta = angle.theta() ;
    EXPECT_NEAR( R1.cos() , cos , tol ) ;                 // cos
    EXPECT_NEAR( R1.sin() , sin , tol ) ;                 // sin
    EXPECT_NEAR( R1.theta() , 2*theta , tol ) ;           // theta
    EXPECT_NEAR( R2.cos() , std::cos( theta2 ) , tol ) ;  // cos
    EXPECT_NEAR( R2.sin() , std::sin( theta2 ) , tol ) ;  // sin
    EXPECT_NEAR( R2.theta() , 2*theta2 , tol ) ;          // theta
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >        angle1( theta1 ) ;
    const qclab::QRotation< R >     rot1( angle1 ) ;
    qclab::qgates::RotationZZ< T >  R1( 0 , 1 , rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >        angle2( theta2 ) ;
    const qclab::QRotation< R >     rot2( angle2 ) ;
    qclab::qgates::RotationZZ< T >  R2( 0 , 1 , rot2 ) ;

    // operator /=
    R1 /= R2 ;
    qclab::QAngle< R >  angle = angle1 - angle2 ;
    R cos   = angle.cos() ;
    R sin   = angle.sin() ;
    R theta = angle.theta() ;
    EXPECT_NEAR( R1.cos() , cos , tol ) ;                 // cos
    EXPECT_NEAR( R1.sin() , sin , tol ) ;                 // sin
    EXPECT_NEAR( R1.theta() , 2*theta , tol ) ;           // theta
    EXPECT_NEAR( R2.cos() , std::cos( theta2 ) , tol ) ;  // cos
    EXPECT_NEAR( R2.sin() , std::sin( theta2 ) , tol ) ;  // sin
    EXPECT_NEAR( R2.theta() , 2*theta2 , tol ) ;          // theta
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >        angle1( theta1 ) ;
    const qclab::QRotation< R >     rot1( angle1 ) ;
    qclab::qgates::RotationZZ< T >  R1( 0 , 1 , rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >        angle2( theta2 ) ;
    const qclab::QRotation< R >     rot2( angle2 ) ;
    qclab::qgates::RotationZZ< T >  R2( 0 , 1 , rot2 ) ;

    // operator *
    qclab::qgates::RotationZZ< T > R12 = R1 * R2 ;
    qclab::QAngle< R >  angle = angle1 + angle2 ;
    R cos   = angle.cos() ;
    R sin   = angle.sin() ;
    R theta = angle.theta() ;
    EXPECT_NEAR( R12.cos() , cos , tol ) ;        // cos
    EXPECT_NEAR( R12.sin() , sin , tol ) ;        // sin
    EXPECT_NEAR( R12.theta() , 2*theta , tol ) ;  // theta
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >        angle1( theta1 ) ;
    const qclab::QRotation< R >     rot1( angle1 ) ;
    qclab::qgates::RotationZZ< T >  R1( 0 , 1 , rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >        angle2( theta2 ) ;
    const qclab::QRotation< R >     rot2( angle2 ) ;
    qclab::qgates::RotationZZ< T >  R2( 0 , 1 , rot2 ) ;

    // operator /
    qclab::qgates::RotationZZ< T > R12 = R1 / R2 ;
    qclab::QAngle< R >  angle = angle1 - angle2 ;
    R cos   = angle.cos() ;
    R sin   = angle.sin() ;
    R theta = angle.theta() ;
    EXPECT_NEAR( R12.cos() , cos , tol ) ;        // cos
    EXPECT_NEAR( R12.sin() , sin , tol ) ;        // sin
    EXPECT_NEAR( R12.theta() , 2*theta , tol ) ;  // theta
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >        angle1( theta1 ) ;
    const qclab::QRotation< R >     rot1( angle1 ) ;
    qclab::qgates::RotationZZ< T >  R1( 0 , 1 , rot1 ) ;
    R cos   = R1.cos() ;
    R sin   = R1.sin() ;
    R theta = R1.theta() ;

    // inv
    qclab::qgates::RotationZZ< T >  R2 = R1.inv() ;
    EXPECT_NEAR( R1.cos() ,  cos , tol ) ;      // cos
    EXPECT_NEAR( R2.cos() ,  cos , tol ) ;      // cos
    EXPECT_NEAR( R1.sin() ,  sin , tol ) ;      // sin
    EXPECT_NEAR( R2.sin() , -sin , tol ) ;      // sin
    EXPECT_NEAR( R1.theta() ,  theta , tol ) ;  // theta
    EXPECT_NEAR( R2.theta() , -theta , tol ) ;  // theta
  }

}


/*
 * complex float
 */
TEST( qclab_qgates_RotationZZ , complex_float ) {
  test_qclab_qgates_RotationZZ< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_RotationZZ , complex_double ) {
  test_qclab_qgates_RotationZZ< std::complex< double > >() ;
}

