#include <gtest/gtest.h>
#include "qgates/CPhase.hpp"

template <typename T>
void test_qclab_qgates_CPhase() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  //
  // CPhases |1> controlled
  //
  {
    qclab::qgates::CPhase< T >  cphase ;

    EXPECT_EQ( cphase.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( cphase.fixed() ) ;          // fixed
    EXPECT_TRUE( cphase.controlled() ) ;      // controlled
    EXPECT_EQ( cphase.control() , 0 ) ;       // control
    EXPECT_EQ( cphase.target() , 1 ) ;        // target
    EXPECT_EQ( cphase.controlState() , 1 ) ;  // controlState

    // matrix
    auto eye = qclab::dense::eye< T >( 4 ) ;
    EXPECT_TRUE( cphase.matrix() == eye ) ;

    // qubit
    EXPECT_EQ( cphase.qubit() , 0 ) ;

    // qubits
    auto qubits = cphase.qubits() ;
    EXPECT_EQ( qubits.size() , 2 ) ;
    EXPECT_EQ( qubits[0] , 0 ) ;
    EXPECT_EQ( qubits[1] , 1 ) ;
    int qnew[] = { 5 , 3 } ;
    cphase.setQubits( &qnew[0] ) ;
    EXPECT_EQ( cphase.qubits()[0] , 3 ) ;
    EXPECT_EQ( cphase.qubits()[1] , 5 ) ;
    qnew[0] = 0 ;
    qnew[1] = 1 ;
    cphase.setQubits( &qnew[0] ) ;

    // print
    cphase.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cphase.toQASM( qasm ) , 0 ) ;
    auto str_theta = qclab::qasm( cphase.theta() ) ;
    std::string qasm_check = "cp(" + str_theta + ") q[0], q[1];\n" ;
    EXPECT_EQ( qasm.str() , qasm_check ) ;
    std::cout << qasm.str() ;

    // gate
    qclab::qgates::Phase< T >  P ;
    EXPECT_TRUE( *cphase.gate() == P ) ;
    EXPECT_TRUE( cphase.gate()->matrix() == P.matrix() ) ;

    // operators == and !=
    EXPECT_TRUE(  cphase != P ) ;
    EXPECT_FALSE( cphase == P ) ;
    qclab::qgates::CPhase< T >  cphase2 ;
    EXPECT_TRUE(  cphase == cphase2 ) ;
    EXPECT_FALSE( cphase != cphase2 ) ;

    // setControl, setTarget, setControlState
    cphase.setControl( 3 ) ;
    EXPECT_EQ( cphase.control() , 3 ) ;
    cphase.setTarget( 5 ) ;
    EXPECT_EQ( cphase.target() , 5 ) ;
    EXPECT_TRUE(  cphase == cphase2 ) ;
    EXPECT_FALSE( cphase != cphase2 ) ;

    cphase.setControl( 4 ) ;
    EXPECT_EQ( cphase.control() , 4 ) ;
    cphase.setTarget( 1 ) ;
    EXPECT_EQ( cphase.target() , 1 ) ;
    EXPECT_TRUE(  cphase != cphase2 ) ;
    EXPECT_FALSE( cphase == cphase2 ) ;

    qubits[0] = 1 ;
    qubits[1] = 2 ;
    cphase.setQubits( &qubits[0] ) ;
    EXPECT_EQ( cphase.control() , 1 ) ;
    EXPECT_EQ( cphase.target() , 2 ) ;
    EXPECT_TRUE(  cphase == cphase2 ) ;
    EXPECT_FALSE( cphase != cphase2 ) ;

    cphase.setControl( 0 ) ;
    EXPECT_EQ( cphase.control() , 0 ) ;
    cphase.setTarget( 1 ) ;
    EXPECT_EQ( cphase.target() , 1 ) ;
    cphase.setControlState( 0 ) ;
    EXPECT_EQ( cphase.controlState() , 0 ) ;
    EXPECT_TRUE(  cphase != cphase2 ) ;
    EXPECT_FALSE( cphase == cphase2 ) ;

    // makeFixed, makeVariable
    cphase.makeFixed() ;
    EXPECT_TRUE( cphase.fixed() ) ;
    cphase.makeVariable() ;
    EXPECT_FALSE( cphase.fixed() ) ;

    // angle, theta, sin, cos
    qclab::QAngle< R > angle ;
    EXPECT_TRUE( cphase.angle() == angle ) ;
    EXPECT_EQ( cphase.theta() , 0 ) ;
    EXPECT_EQ( cphase.cos() , 1 ) ;
    EXPECT_EQ( cphase.sin() , 0 ) ;

    // update(angle)
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    angle.update( pi/4 ) ;
    cphase.update( angle ) ;
    EXPECT_NEAR( cphase.theta() , pi/4 , eps ) ;
    EXPECT_NEAR( cphase.cos() , cos , eps ) ;
    EXPECT_NEAR( cphase.sin() , sin , eps ) ;

    // update(theta)
    cphase.update( pi/2 ) ;
    EXPECT_NEAR( cphase.theta() , pi/2 , eps ) ;
    EXPECT_NEAR( cphase.cos() , 0 , eps ) ;
    EXPECT_NEAR( cphase.sin() , 1 , eps ) ;

    // update(cos,sin)
    cphase.update( cos , sin ) ;
    EXPECT_NEAR( cphase.theta() , pi/4 , eps ) ;
    EXPECT_NEAR( cphase.cos() , cos , eps ) ;
    EXPECT_NEAR( cphase.sin() , sin , eps ) ;
  }

  {
    const R theta = pi/4 ;
    qclab::QAngle< R >  angle( theta ) ;
    qclab::qgates::CPhase< T >  cphase( 1 , 3 , angle ) ;

    EXPECT_EQ( cphase.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( cphase.fixed() ) ;          // fixed
    EXPECT_TRUE( cphase.controlled() ) ;      // controlled
    EXPECT_EQ( cphase.control() , 1 ) ;       // control
    EXPECT_EQ( cphase.target() , 3 ) ;        // target
    EXPECT_EQ( cphase.controlState() , 1 ) ;  // controlState

    EXPECT_NEAR( cphase.theta() , theta , eps ) ;
    EXPECT_NEAR( cphase.cos() , std::cos( theta ) , eps ) ;
    EXPECT_NEAR( cphase.sin() , std::sin( theta ) , eps ) ;
  }

  {
    const R theta = pi/4 ;
    qclab::qgates::CPhase< T >  cphase( 1 , 3 , theta ) ;

    EXPECT_EQ( cphase.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( cphase.fixed() ) ;          // fixed
    EXPECT_TRUE( cphase.controlled() ) ;      // controlled
    EXPECT_EQ( cphase.control() , 1 ) ;       // control
    EXPECT_EQ( cphase.target() , 3 ) ;        // target
    EXPECT_EQ( cphase.controlState() , 1 ) ;  // controlState

    EXPECT_NEAR( cphase.theta() , theta , eps ) ;
    EXPECT_NEAR( cphase.cos() , std::cos( theta ) , eps ) ;
    EXPECT_NEAR( cphase.sin() , std::sin( theta ) , eps ) ;
  }

  {
    const R theta = pi/4 ;
    const R cos = std::cos( theta ) ;
    const R sin = std::sin( theta ) ;
    qclab::qgates::CPhase< T >  cphase( 1 , 3 , cos , sin ) ;

    EXPECT_EQ( cphase.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( cphase.fixed() ) ;          // fixed
    EXPECT_TRUE( cphase.controlled() ) ;      // controlled
    EXPECT_EQ( cphase.control() , 1 ) ;       // control
    EXPECT_EQ( cphase.target() , 3 ) ;        // target
    EXPECT_EQ( cphase.controlState() , 1 ) ;  // controlState

    EXPECT_NEAR( cphase.theta() , theta , eps ) ;
    EXPECT_NEAR( cphase.cos() , cos , eps ) ;
    EXPECT_NEAR( cphase.sin() , sin , eps ) ;
  }

}


/*
 * complex float
 */
TEST( qclab_qgates_CPhase , complex_float ) {
  test_qclab_qgates_CPhase< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_CPhase , complex_double ) {
  test_qclab_qgates_CPhase< std::complex< double > >() ;
}

