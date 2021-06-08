#include <gtest/gtest.h>
#include "qgates/CRotationX.hpp"

template <typename T>
void test_qclab_qgates_CRotationX() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  //
  // CRotationXs |1> controlled
  //
  {
    qclab::qgates::CRotationX< T >  crotx ;

    EXPECT_EQ( crotx.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( crotx.fixed() ) ;          // fixed
    EXPECT_TRUE( crotx.controlled() ) ;      // controlled
    EXPECT_EQ( crotx.control() , 0 ) ;       // control
    EXPECT_EQ( crotx.target() , 1 ) ;        // target
    EXPECT_EQ( crotx.controlState() , 1 ) ;  // controlState

    // matrix
    auto eye = qclab::dense::eye< T >( 4 ) ;
    EXPECT_TRUE( crotx.matrix() == eye ) ;

    // qubit
    EXPECT_EQ( crotx.qubit() , 0 ) ;

    // qubits
    auto qubits = crotx.qubits() ;
    EXPECT_EQ( qubits.size() , 2 ) ;
    EXPECT_EQ( qubits[0] , 0 ) ;
    EXPECT_EQ( qubits[1] , 1 ) ;
    int qnew[] = { 5 , 3 } ;
    crotx.setQubits( &qnew[0] ) ;
    EXPECT_EQ( crotx.qubits()[0] , 3 ) ;
    EXPECT_EQ( crotx.qubits()[1] , 5 ) ;
    qnew[0] = 0 ;
    qnew[1] = 1 ;
    crotx.setQubits( &qnew[0] ) ;

    // print
    crotx.update( pi/2 ) ;
    crotx.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( crotx.toQASM( qasm ) , 0 ) ;
    auto str_theta = qclab::qasm( crotx.theta() ) ;
    std::string qasm_check = "crx(" + str_theta + ") q[0], q[1];\n" ;
    EXPECT_EQ( qasm.str() , qasm_check ) ;
    std::cout << qasm.str() ;
    crotx.update( 0. ) ;

    // gate
    qclab::qgates::RotationX< T >  P ;
    EXPECT_TRUE( *crotx.gate() == P ) ;
    EXPECT_TRUE( crotx.gate()->matrix() == P.matrix() ) ;

    // operators == and !=
    EXPECT_TRUE(  crotx != P ) ;
    EXPECT_FALSE( crotx == P ) ;
    qclab::qgates::CRotationX< T >  crotx2 ;
    EXPECT_TRUE(  crotx == crotx2 ) ;
    EXPECT_FALSE( crotx != crotx2 ) ;

    // setControl, setTarget, setControlState
    crotx.setControl( 3 ) ;
    EXPECT_EQ( crotx.control() , 3 ) ;
    crotx.setTarget( 5 ) ;
    EXPECT_EQ( crotx.target() , 5 ) ;
    EXPECT_TRUE(  crotx == crotx2 ) ;
    EXPECT_FALSE( crotx != crotx2 ) ;

    crotx.setControl( 4 ) ;
    EXPECT_EQ( crotx.control() , 4 ) ;
    crotx.setTarget( 1 ) ;
    EXPECT_EQ( crotx.target() , 1 ) ;
    EXPECT_TRUE(  crotx != crotx2 ) ;
    EXPECT_FALSE( crotx == crotx2 ) ;

    qubits[0] = 1 ;
    qubits[1] = 2 ;
    crotx.setQubits( &qubits[0] ) ;
    EXPECT_EQ( crotx.control() , 1 ) ;
    EXPECT_EQ( crotx.target() , 2 ) ;
    EXPECT_TRUE(  crotx == crotx2 ) ;
    EXPECT_FALSE( crotx != crotx2 ) ;

    crotx.setControl( 0 ) ;
    EXPECT_EQ( crotx.control() , 0 ) ;
    crotx.setTarget( 1 ) ;
    EXPECT_EQ( crotx.target() , 1 ) ;
    crotx.setControlState( 0 ) ;
    EXPECT_EQ( crotx.controlState() , 0 ) ;
    EXPECT_TRUE(  crotx != crotx2 ) ;
    EXPECT_FALSE( crotx == crotx2 ) ;

    // makeFixed, makeVariable
    crotx.makeFixed() ;
    EXPECT_TRUE( crotx.fixed() ) ;
    crotx.makeVariable() ;
    EXPECT_FALSE( crotx.fixed() ) ;

    // angle, theta, sin, cos
    qclab::QAngle< R > angle ;
    EXPECT_TRUE( crotx.angle() == angle ) ;
    EXPECT_EQ( crotx.theta() , 0 ) ;
    EXPECT_EQ( crotx.cos() , 1 ) ;
    EXPECT_EQ( crotx.sin() , 0 ) ;

    // update(angle)
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    angle.update( pi/4 ) ;
    crotx.update( angle ) ;
    EXPECT_NEAR( crotx.theta() , pi/2 , eps ) ;
    EXPECT_NEAR( crotx.cos() , cos , eps ) ;
    EXPECT_NEAR( crotx.sin() , sin , eps ) ;

    // update(theta)
    crotx.update( pi ) ;
    EXPECT_NEAR( crotx.theta() , pi , eps ) ;
    EXPECT_NEAR( crotx.cos() , 0 , eps ) ;
    EXPECT_NEAR( crotx.sin() , 1 , eps ) ;

    // update(cos,sin)
    crotx.update( cos , sin ) ;
    EXPECT_NEAR( crotx.theta() , pi/2 , eps ) ;
    EXPECT_NEAR( crotx.cos() , cos , eps ) ;
    EXPECT_NEAR( crotx.sin() , sin , eps ) ;
  }

  {
    const R theta = pi/4 ;
    qclab::QAngle< R >  angle( theta/2 ) ;
    qclab::qgates::CRotationX< T >  crotx( 1 , 3 , angle ) ;

    EXPECT_EQ( crotx.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( crotx.fixed() ) ;          // fixed
    EXPECT_TRUE( crotx.controlled() ) ;      // controlled
    EXPECT_EQ( crotx.control() , 1 ) ;       // control
    EXPECT_EQ( crotx.target() , 3 ) ;        // target
    EXPECT_EQ( crotx.controlState() , 1 ) ;  // controlState

    EXPECT_NEAR( crotx.theta() , theta , eps ) ;
    EXPECT_NEAR( crotx.cos() , std::cos( theta/2 ) , eps ) ;
    EXPECT_NEAR( crotx.sin() , std::sin( theta/2 ) , eps ) ;
  }

  {
    const R theta = pi/4 ;
    qclab::qgates::CRotationX< T >  crotx( 1 , 3 , theta ) ;

    EXPECT_EQ( crotx.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( crotx.fixed() ) ;          // fixed
    EXPECT_TRUE( crotx.controlled() ) ;      // controlled
    EXPECT_EQ( crotx.control() , 1 ) ;       // control
    EXPECT_EQ( crotx.target() , 3 ) ;        // target
    EXPECT_EQ( crotx.controlState() , 1 ) ;  // controlState

    EXPECT_NEAR( crotx.theta() , theta , eps ) ;
    EXPECT_NEAR( crotx.cos() , std::cos( theta/2 ) , eps ) ;
    EXPECT_NEAR( crotx.sin() , std::sin( theta/2 ) , eps ) ;
  }

  {
    const R theta = pi/4 ;
    const R cos = std::cos( theta ) ;
    const R sin = std::sin( theta ) ;
    qclab::qgates::CRotationX< T >  crotx( 1 , 3 , cos , sin ) ;

    EXPECT_EQ( crotx.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( crotx.fixed() ) ;          // fixed
    EXPECT_TRUE( crotx.controlled() ) ;      // controlled
    EXPECT_EQ( crotx.control() , 1 ) ;       // control
    EXPECT_EQ( crotx.target() , 3 ) ;        // target
    EXPECT_EQ( crotx.controlState() , 1 ) ;  // controlState

    EXPECT_NEAR( crotx.theta() , 2*theta , eps ) ;
    EXPECT_NEAR( crotx.cos() , cos , eps ) ;
    EXPECT_NEAR( crotx.sin() , sin , eps ) ;
  }

}


/*
 * complex float
 */
TEST( qclab_qgates_CRotationX , complex_float ) {
  test_qclab_qgates_CRotationX< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_CRotationX , complex_double ) {
  test_qclab_qgates_CRotationX< std::complex< double > >() ;
}

