#include <gtest/gtest.h>
#include "qclab/qgates/RotationXX.hpp"

template <typename T>
void test_qclab_qgates_RotationXX() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R tol = 10 * std::numeric_limits< R >::epsilon() ;

  {
    qclab::qgates::RotationXX< T >  Rxx ;

    EXPECT_EQ( Rxx.nbQubits() , 2 ) ;   // nbQubits
    EXPECT_FALSE( Rxx.fixed() ) ;       // fixed
    EXPECT_FALSE( Rxx.controlled() ) ;  // controlled
    EXPECT_EQ( Rxx.cos() , 1.0 ) ;      // cos
    EXPECT_EQ( Rxx.sin() , 0.0 ) ;      // sin
    EXPECT_EQ( Rxx.theta() , 0.0 ) ;    // theta

    // matrix
    auto eye = qclab::dense::eye< T >( 4 ) ;
    EXPECT_TRUE( Rxx.matrix() == eye ) ;

    // qubits
    EXPECT_EQ( Rxx.qubits().size() , 2 ) ;
    EXPECT_EQ( Rxx.qubits()[0] , 0 ) ;
    EXPECT_EQ( Rxx.qubits()[1] , 1 ) ;
    int qnew[2] = { 3 , 5 } ;
    Rxx.setQubits( &qnew[0] ) ;
    EXPECT_EQ( Rxx.qubits().size() , 2 ) ;
    EXPECT_EQ( Rxx.qubits()[0] , 3 ) ;
    EXPECT_EQ( Rxx.qubits()[1] , 5 ) ;

    // fixed
    Rxx.makeFixed() ;
    EXPECT_TRUE( Rxx.fixed() ) ;
    Rxx.makeVariable() ;
    EXPECT_FALSE( Rxx.fixed() ) ;

    // update(rot)
    qclab::QRotation< R >  new_rot( 1 ) ;
    Rxx.update( new_rot ) ;
    EXPECT_NEAR( Rxx.theta() , 1 , tol ) ;
    EXPECT_NEAR( Rxx.cos() , std::cos( 0.5 ) , tol ) ;
    EXPECT_NEAR( Rxx.sin() , std::sin( 0.5 ) , tol ) ;

    // update(theta)
    Rxx.update( pi/2 ) ;
    EXPECT_NEAR( Rxx.theta() , pi/2 , tol ) ;
    EXPECT_NEAR( Rxx.cos() , std::cos( pi/4 ) , tol ) ;
    EXPECT_NEAR( Rxx.sin() , std::sin( pi/4 ) , tol ) ;

    // matrix
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    EXPECT_EQ( Rxx.matrix()(0,0) , T(cos,0.0) ) ;
    EXPECT_EQ( Rxx.matrix()(1,1) , T(cos,0.0) ) ;
    EXPECT_EQ( Rxx.matrix()(2,2) , T(cos,0.0) ) ;
    EXPECT_EQ( Rxx.matrix()(3,3) , T(cos,0.0) ) ;
    EXPECT_EQ( Rxx.matrix()(3,0) , T(0.0,-sin) ) ;
    EXPECT_EQ( Rxx.matrix()(2,1) , T(0.0,-sin) ) ;
    EXPECT_EQ( Rxx.matrix()(1,2) , T(0.0,-sin) ) ;
    EXPECT_EQ( Rxx.matrix()(0,3) , T(0.0,-sin) ) ;
    EXPECT_EQ( Rxx.matrix()(1,0) , T(0) ) ;
    EXPECT_EQ( Rxx.matrix()(2,0) , T(0) ) ;
    EXPECT_EQ( Rxx.matrix()(0,1) , T(0) ) ;
    EXPECT_EQ( Rxx.matrix()(3,1) , T(0) ) ;
    EXPECT_EQ( Rxx.matrix()(0,2) , T(0) ) ;
    EXPECT_EQ( Rxx.matrix()(3,2) , T(0) ) ;
    EXPECT_EQ( Rxx.matrix()(1,3) , T(0) ) ;
    EXPECT_EQ( Rxx.matrix()(2,3) , T(0) ) ;

    // print
    Rxx.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( Rxx.toQASM( qasm ) , 0 ) ;
    auto str_theta = qclab::qasm( Rxx.theta() ) ;
    std::string qasm_check = "rxx(" + str_theta + ") q[3], q[5];\n" ;
    EXPECT_EQ( qasm.str() , qasm_check ) ;
    std::cout << qasm.str() ;

    // update(cos,sin)
    Rxx.update( std::cos( pi/3 ) , std::sin( pi/3 ) ) ;
    EXPECT_NEAR( Rxx.cos() , std::cos( pi/3 ) , tol ) ;
    EXPECT_NEAR( Rxx.sin() , std::sin( pi/3 ) , tol ) ;
    EXPECT_NEAR( Rxx.theta() , 2*(pi/3) , tol ) ;

    // operators == and !=
    qclab::qgates::RotationXX< T >  Rxx2( std::cos( pi/3 ) , std::sin( pi/3 ) );
    EXPECT_TRUE( Rxx == Rxx2 ) ;
    EXPECT_FALSE( Rxx != Rxx2 ) ;
    Rxx2.update( 1 ) ;
    EXPECT_TRUE( Rxx != Rxx2 ) ;
    EXPECT_FALSE( Rxx == Rxx2 ) ;
  }

  //
  // constructors without qubits
  //
  {
    qclab::QRotation< R >           rot( pi/2 ) ;
    qclab::qgates::RotationXX< T >  Rxx( rot ) ;

    EXPECT_EQ( Rxx.nbQubits() , 2 ) ;                    // nbQubits
    EXPECT_FALSE( Rxx.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rxx.controlled() ) ;                   // controlled

    auto qubits = Rxx.qubits() ;
    EXPECT_EQ( qubits[0] , 0 ) ;                         // qubit0
    EXPECT_EQ( qubits[1] , 1 ) ;                         // qubit1

    EXPECT_NEAR( Rxx.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Rxx.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Rxx.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    qclab::qgates::RotationXX< T >  Rxx( pi/2 ) ;

    EXPECT_EQ( Rxx.nbQubits() , 2 ) ;                    // nbQubits
    EXPECT_FALSE( Rxx.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rxx.controlled() ) ;                   // controlled

    auto qubits = Rxx.qubits() ;
    EXPECT_EQ( qubits[0] , 0 ) ;                         // qubit0
    EXPECT_EQ( qubits[1] , 1 ) ;                         // qubit1

    EXPECT_NEAR( Rxx.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Rxx.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Rxx.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationXX< T >  Rxx( cos , sin ) ;

    EXPECT_EQ( Rxx.nbQubits() , 2 ) ;          // nbQubits
    EXPECT_FALSE( Rxx.fixed() ) ;              // fixed
    EXPECT_FALSE( Rxx.controlled() ) ;         // controlled

    auto qubits = Rxx.qubits() ;
    EXPECT_EQ( qubits[0] , 0 ) ;               // qubit0
    EXPECT_EQ( qubits[1] , 1 ) ;               // qubit1

    EXPECT_NEAR( Rxx.theta() , pi/2 , tol ) ;  // theta
    EXPECT_NEAR( Rxx.cos() , cos , tol ) ;     // cos
    EXPECT_NEAR( Rxx.sin() , sin , tol ) ;     // sin
  }


  //
  // constructors with *qubits
  //
  {
    int qbits[2] = { 3 , 5 } ;
    qclab::QRotation< R >           rot( pi/2 ) ;
    qclab::qgates::RotationXX< T >  Rxx( &qbits[0] , rot ) ;

    EXPECT_EQ( Rxx.nbQubits() , 2 ) ;                    // nbQubits
    EXPECT_FALSE( Rxx.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rxx.controlled() ) ;                   // controlled

    auto qubits = Rxx.qubits() ;
    EXPECT_EQ( qubits[0] , 3 ) ;                         // qubit0
    EXPECT_EQ( qubits[1] , 5 ) ;                         // qubit1

    EXPECT_NEAR( Rxx.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Rxx.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Rxx.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    int qbits[2] = { 3 , 5 } ;
    qclab::qgates::RotationXX< T >  Rxx( &qbits[0] , pi/2 ) ;

    EXPECT_EQ( Rxx.nbQubits() , 2 ) ;                    // nbQubits
    EXPECT_FALSE( Rxx.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rxx.controlled() ) ;                   // controlled

    auto qubits = Rxx.qubits() ;
    EXPECT_EQ( qubits[0] , 3 ) ;                         // qubit0
    EXPECT_EQ( qubits[1] , 5 ) ;                         // qubit1

    EXPECT_NEAR( Rxx.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Rxx.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Rxx.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    int qbits[2] = { 3 , 5 } ;
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationXX< T >  Rxx( &qbits[0] , cos , sin ) ;

    EXPECT_EQ( Rxx.nbQubits() , 2 ) ;          // nbQubits
    EXPECT_FALSE( Rxx.fixed() ) ;              // fixed
    EXPECT_FALSE( Rxx.controlled() ) ;         // controlled

    auto qubits = Rxx.qubits() ;
    EXPECT_EQ( qubits[0] , 3 ) ;               // qubit0
    EXPECT_EQ( qubits[1] , 5 ) ;               // qubit1

    EXPECT_NEAR( Rxx.theta() , pi/2 , tol ) ;  // theta
    EXPECT_NEAR( Rxx.cos() , cos , tol ) ;     // cos
    EXPECT_NEAR( Rxx.sin() , sin , tol ) ;     // sin
  }


  //
  // constructors with qubit0 and qubit1
  //
  {
    int qubit0 = 3 ;
    int qubit1 = 5 ;
    qclab::QRotation< R >           rot( pi/2 ) ;
    qclab::qgates::RotationXX< T >  Rxx( qubit0 , qubit1 , rot ) ;

    EXPECT_EQ( Rxx.nbQubits() , 2 ) ;                    // nbQubits
    EXPECT_FALSE( Rxx.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rxx.controlled() ) ;                   // controlled

    auto qubits = Rxx.qubits() ;
    EXPECT_EQ( qubits[0] , 3 ) ;                         // qubit0
    EXPECT_EQ( qubits[1] , 5 ) ;                         // qubit1

    EXPECT_NEAR( Rxx.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Rxx.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Rxx.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    int qubit0 = 3 ;
    int qubit1 = 5 ;
    qclab::qgates::RotationXX< T >  Rxx( qubit0 , qubit1 , pi/2 ) ;

    EXPECT_EQ( Rxx.nbQubits() , 2 ) ;                    // nbQubits
    EXPECT_FALSE( Rxx.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rxx.controlled() ) ;                   // controlled

    auto qubits = Rxx.qubits() ;
    EXPECT_EQ( qubits[0] , 3 ) ;                         // qubit0
    EXPECT_EQ( qubits[1] , 5 ) ;                         // qubit1

    EXPECT_NEAR( Rxx.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Rxx.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Rxx.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    int qubit0 = 3 ;
    int qubit1 = 5 ;
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationXX< T >  Rxx( qubit0 , qubit1 , cos , sin ) ;

    EXPECT_EQ( Rxx.nbQubits() , 2 ) ;          // nbQubits
    EXPECT_FALSE( Rxx.fixed() ) ;              // fixed
    EXPECT_FALSE( Rxx.controlled() ) ;         // controlled

    auto qubits = Rxx.qubits() ;
    EXPECT_EQ( qubits[0] , 3 ) ;               // qubit0
    EXPECT_EQ( qubits[1] , 5 ) ;               // qubit1

    EXPECT_NEAR( Rxx.theta() , pi/2 , tol ) ;  // theta
    EXPECT_NEAR( Rxx.cos() , cos , tol ) ;     // cos
    EXPECT_NEAR( Rxx.sin() , sin , tol ) ;     // sin
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >        angle1( theta1 ) ;
    const qclab::QRotation< R >     rot1( angle1 ) ;
    qclab::qgates::RotationXX< T >  R1( 0 , 1 , rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >        angle2( theta2 ) ;
    const qclab::QRotation< R >     rot2( angle2 ) ;
    qclab::qgates::RotationXX< T >  R2( 0 , 1 , rot2 ) ;

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
    qclab::qgates::RotationXX< T >  R1( 0 , 1 , rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >        angle2( theta2 ) ;
    const qclab::QRotation< R >     rot2( angle2 ) ;
    qclab::qgates::RotationXX< T >  R2( 0 , 1 , rot2 ) ;

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
    qclab::qgates::RotationXX< T >  R1( 0 , 1 , rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >        angle2( theta2 ) ;
    const qclab::QRotation< R >     rot2( angle2 ) ;
    qclab::qgates::RotationXX< T >  R2( 0 , 1 , rot2 ) ;

    // operator *
    qclab::qgates::RotationXX< T > R12 = R1 * R2 ;
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
    qclab::qgates::RotationXX< T >  R1( 0 , 1 , rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >        angle2( theta2 ) ;
    const qclab::QRotation< R >     rot2( angle2 ) ;
    qclab::qgates::RotationXX< T >  R2( 0 , 1 , rot2 ) ;

    // operator /
    qclab::qgates::RotationXX< T > R12 = R1 / R2 ;
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
    qclab::qgates::RotationXX< T >  R1( 0 , 1 , rot1 ) ;
    R cos   = R1.cos() ;
    R sin   = R1.sin() ;
    R theta = R1.theta() ;

    // inv
    qclab::qgates::RotationXX< T >  R2 = R1.inv() ;
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
TEST( qclab_qgates_RotationXX , complex_float ) {
  test_qclab_qgates_RotationXX< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_RotationXX , complex_double ) {
  test_qclab_qgates_RotationXX< std::complex< double > >() ;
}

