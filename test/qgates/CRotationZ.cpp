#include <gtest/gtest.h>
#include "qclab/qgates/CRotationZ.hpp"

template <typename T>
void test_qclab_qgates_CRotationZ() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  //
  // CRotationZs |1> controlled
  //
  {
    qclab::qgates::CRotationZ< T >  crotz ;

    EXPECT_EQ( crotz.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( crotz.fixed() ) ;          // fixed
    EXPECT_TRUE( crotz.controlled() ) ;      // controlled
    EXPECT_EQ( crotz.control() , 0 ) ;       // control
    EXPECT_EQ( crotz.target() , 1 ) ;        // target
    EXPECT_EQ( crotz.controlState() , 1 ) ;  // controlState

    // matrix
    auto eye = qclab::dense::eye< T >( 4 ) ;
    EXPECT_TRUE( crotz.matrix() == eye ) ;

    // qubit
    EXPECT_EQ( crotz.qubit() , 0 ) ;

    // qubits
    auto qubits = crotz.qubits() ;
    EXPECT_EQ( qubits.size() , 2 ) ;
    EXPECT_EQ( qubits[0] , 0 ) ;
    EXPECT_EQ( qubits[1] , 1 ) ;
    int qnew[] = { 5 , 3 } ;
    crotz.setQubits( &qnew[0] ) ;
    EXPECT_EQ( crotz.qubits()[0] , 3 ) ;
    EXPECT_EQ( crotz.qubits()[1] , 5 ) ;
    qnew[0] = 0 ;
    qnew[1] = 1 ;
    crotz.setQubits( &qnew[0] ) ;

    // print
    crotz.update( pi/2 ) ;
    crotz.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( crotz.toQASM( qasm ) , 0 ) ;
    auto str_theta = qclab::qasm( crotz.theta() ) ;
    std::string qasm_check = "crz(" + str_theta + ") q[0], q[1];\n" ;
    EXPECT_EQ( qasm.str() , qasm_check ) ;
    std::cout << qasm.str() ;
    crotz.update( 0. ) ;

    // gate
    qclab::qgates::RotationZ< T >  P ;
    EXPECT_TRUE( *crotz.gate() == P ) ;
    EXPECT_TRUE( crotz.gate()->matrix() == P.matrix() ) ;

    // operators == and !=
    EXPECT_TRUE(  crotz != P ) ;
    EXPECT_FALSE( crotz == P ) ;
    qclab::qgates::CRotationZ< T >  crotz2 ;
    EXPECT_TRUE(  crotz == crotz2 ) ;
    EXPECT_FALSE( crotz != crotz2 ) ;

    // setControl, setTarget, setControlState
    crotz.setControl( 3 ) ;
    EXPECT_EQ( crotz.control() , 3 ) ;
    crotz.setTarget( 5 ) ;
    EXPECT_EQ( crotz.target() , 5 ) ;
    EXPECT_TRUE(  crotz == crotz2 ) ;
    EXPECT_FALSE( crotz != crotz2 ) ;

    crotz.setControl( 4 ) ;
    EXPECT_EQ( crotz.control() , 4 ) ;
    crotz.setTarget( 1 ) ;
    EXPECT_EQ( crotz.target() , 1 ) ;
    EXPECT_TRUE(  crotz != crotz2 ) ;
    EXPECT_FALSE( crotz == crotz2 ) ;

    qubits[0] = 1 ;
    qubits[1] = 2 ;
    crotz.setQubits( &qubits[0] ) ;
    EXPECT_EQ( crotz.control() , 1 ) ;
    EXPECT_EQ( crotz.target() , 2 ) ;
    EXPECT_TRUE(  crotz == crotz2 ) ;
    EXPECT_FALSE( crotz != crotz2 ) ;

    crotz.setControl( 0 ) ;
    EXPECT_EQ( crotz.control() , 0 ) ;
    crotz.setTarget( 1 ) ;
    EXPECT_EQ( crotz.target() , 1 ) ;
    crotz.setControlState( 0 ) ;
    EXPECT_EQ( crotz.controlState() , 0 ) ;
    EXPECT_TRUE(  crotz != crotz2 ) ;
    EXPECT_FALSE( crotz == crotz2 ) ;

    // makeFixed, makeVariable
    crotz.makeFixed() ;
    EXPECT_TRUE( crotz.fixed() ) ;
    crotz.makeVariable() ;
    EXPECT_FALSE( crotz.fixed() ) ;

    // rotation, theta, sin, cos
    qclab::QRotation< R > rot ;
    EXPECT_TRUE( crotz.rotation() == rot ) ;
    EXPECT_EQ( crotz.theta() , 0 ) ;
    EXPECT_EQ( crotz.cos() , 1 ) ;
    EXPECT_EQ( crotz.sin() , 0 ) ;

    // update(rot)
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    rot.update( pi/2 ) ;
    crotz.update( rot ) ;
    EXPECT_NEAR( crotz.theta() , pi/2 , eps ) ;
    EXPECT_NEAR( crotz.cos() , cos , eps ) ;
    EXPECT_NEAR( crotz.sin() , sin , eps ) ;

    // update(theta)
    crotz.update( pi ) ;
    EXPECT_NEAR( crotz.theta() , pi , eps ) ;
    EXPECT_NEAR( crotz.cos() , 0 , eps ) ;
    EXPECT_NEAR( crotz.sin() , 1 , eps ) ;

    // update(cos,sin)
    crotz.update( cos , sin ) ;
    EXPECT_NEAR( crotz.theta() , pi/2 , eps ) ;
    EXPECT_NEAR( crotz.cos() , cos , eps ) ;
    EXPECT_NEAR( crotz.sin() , sin , eps ) ;
  }

  {
    const R theta = pi/4 ;
    qclab::QRotation< R >  rot( theta ) ;
    qclab::qgates::CRotationZ< T >  crotz( 1 , 3 , rot ) ;

    EXPECT_EQ( crotz.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( crotz.fixed() ) ;          // fixed
    EXPECT_TRUE( crotz.controlled() ) ;      // controlled
    EXPECT_EQ( crotz.control() , 1 ) ;       // control
    EXPECT_EQ( crotz.target() , 3 ) ;        // target
    EXPECT_EQ( crotz.controlState() , 1 ) ;  // controlState

    EXPECT_NEAR( crotz.theta() , theta , eps ) ;
    EXPECT_NEAR( crotz.cos() , std::cos( theta/2 ) , eps ) ;
    EXPECT_NEAR( crotz.sin() , std::sin( theta/2 ) , eps ) ;
  }

  {
    const R theta = pi/4 ;
    qclab::qgates::CRotationZ< T >  crotz( 1 , 3 , theta ) ;

    EXPECT_EQ( crotz.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( crotz.fixed() ) ;          // fixed
    EXPECT_TRUE( crotz.controlled() ) ;      // controlled
    EXPECT_EQ( crotz.control() , 1 ) ;       // control
    EXPECT_EQ( crotz.target() , 3 ) ;        // target
    EXPECT_EQ( crotz.controlState() , 1 ) ;  // controlState

    EXPECT_NEAR( crotz.theta() , theta , eps ) ;
    EXPECT_NEAR( crotz.cos() , std::cos( theta/2 ) , eps ) ;
    EXPECT_NEAR( crotz.sin() , std::sin( theta/2 ) , eps ) ;
  }

  {
    const R theta = pi/4 ;
    const R cos = std::cos( theta ) ;
    const R sin = std::sin( theta ) ;
    qclab::qgates::CRotationZ< T >  crotz( 1 , 3 , cos , sin ) ;

    EXPECT_EQ( crotz.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( crotz.fixed() ) ;          // fixed
    EXPECT_TRUE( crotz.controlled() ) ;      // controlled
    EXPECT_EQ( crotz.control() , 1 ) ;       // control
    EXPECT_EQ( crotz.target() , 3 ) ;        // target
    EXPECT_EQ( crotz.controlState() , 1 ) ;  // controlState

    EXPECT_NEAR( crotz.theta() , 2*theta , eps ) ;
    EXPECT_NEAR( crotz.cos() , cos , eps ) ;
    EXPECT_NEAR( crotz.sin() , sin , eps ) ;
  }

}


/*
 * complex float
 */
TEST( qclab_qgates_CRotationZ , complex_float ) {
  test_qclab_qgates_CRotationZ< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_CRotationZ , complex_double ) {
  test_qclab_qgates_CRotationZ< std::complex< double > >() ;
}

