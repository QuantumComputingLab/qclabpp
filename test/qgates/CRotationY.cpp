#include <gtest/gtest.h>
#include "qclab/qgates/CRotationY.hpp"

template <typename T>
void test_qclab_qgates_CRotationY() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  //
  // CRotationYs |1> controlled
  //
  {
    qclab::qgates::CRotationY< T >  croty ;

    EXPECT_EQ( croty.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( croty.fixed() ) ;          // fixed
    EXPECT_TRUE( croty.controlled() ) ;      // controlled
    EXPECT_EQ( croty.control() , 0 ) ;       // control
    EXPECT_EQ( croty.target() , 1 ) ;        // target
    EXPECT_EQ( croty.controlState() , 1 ) ;  // controlState

    // matrix
    auto eye = qclab::dense::eye< T >( 4 ) ;
    EXPECT_TRUE( croty.matrix() == eye ) ;

    // qubit
    EXPECT_EQ( croty.qubit() , 0 ) ;

    // qubits
    auto qubits = croty.qubits() ;
    EXPECT_EQ( qubits.size() , 2 ) ;
    EXPECT_EQ( qubits[0] , 0 ) ;
    EXPECT_EQ( qubits[1] , 1 ) ;
    int qnew[] = { 5 , 3 } ;
    croty.setQubits( &qnew[0] ) ;
    EXPECT_EQ( croty.qubits()[0] , 3 ) ;
    EXPECT_EQ( croty.qubits()[1] , 5 ) ;
    qnew[0] = 0 ;
    qnew[1] = 1 ;
    croty.setQubits( &qnew[0] ) ;

    // print
    croty.update( pi/2 ) ;
    croty.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( croty.toQASM( qasm ) , 0 ) ;
    auto str_theta = qclab::qasm( croty.theta() ) ;
    std::string qasm_check = "cry(" + str_theta + ") q[0], q[1];\n" ;
    EXPECT_EQ( qasm.str() , qasm_check ) ;
    std::cout << qasm.str() ;
    croty.update( 0. ) ;

    // gate
    qclab::qgates::RotationY< T >  P ;
    EXPECT_TRUE( *croty.gate() == P ) ;
    EXPECT_TRUE( croty.gate()->matrix() == P.matrix() ) ;

    // operators == and !=
    EXPECT_TRUE(  croty != P ) ;
    EXPECT_FALSE( croty == P ) ;
    qclab::qgates::CRotationY< T >  croty2 ;
    EXPECT_TRUE(  croty == croty2 ) ;
    EXPECT_FALSE( croty != croty2 ) ;

    // setControl, setTarget, setControlState
    croty.setControl( 3 ) ;
    EXPECT_EQ( croty.control() , 3 ) ;
    croty.setTarget( 5 ) ;
    EXPECT_EQ( croty.target() , 5 ) ;
    EXPECT_TRUE(  croty == croty2 ) ;
    EXPECT_FALSE( croty != croty2 ) ;

    croty.setControl( 4 ) ;
    EXPECT_EQ( croty.control() , 4 ) ;
    croty.setTarget( 1 ) ;
    EXPECT_EQ( croty.target() , 1 ) ;
    EXPECT_TRUE(  croty != croty2 ) ;
    EXPECT_FALSE( croty == croty2 ) ;

    qubits[0] = 1 ;
    qubits[1] = 2 ;
    croty.setQubits( &qubits[0] ) ;
    EXPECT_EQ( croty.control() , 1 ) ;
    EXPECT_EQ( croty.target() , 2 ) ;
    EXPECT_TRUE(  croty == croty2 ) ;
    EXPECT_FALSE( croty != croty2 ) ;

    croty.setControl( 0 ) ;
    EXPECT_EQ( croty.control() , 0 ) ;
    croty.setTarget( 1 ) ;
    EXPECT_EQ( croty.target() , 1 ) ;
    croty.setControlState( 0 ) ;
    EXPECT_EQ( croty.controlState() , 0 ) ;
    EXPECT_TRUE(  croty != croty2 ) ;
    EXPECT_FALSE( croty == croty2 ) ;

    // makeFixed, makeVariable
    croty.makeFixed() ;
    EXPECT_TRUE( croty.fixed() ) ;
    croty.makeVariable() ;
    EXPECT_FALSE( croty.fixed() ) ;

    // angle, theta, sin, cos
    qclab::QAngle< R > angle ;
    EXPECT_TRUE( croty.angle() == angle ) ;
    EXPECT_EQ( croty.theta() , 0 ) ;
    EXPECT_EQ( croty.cos() , 1 ) ;
    EXPECT_EQ( croty.sin() , 0 ) ;

    // update(angle)
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    angle.update( pi/4 ) ;
    croty.update( angle ) ;
    EXPECT_NEAR( croty.theta() , pi/2 , eps ) ;
    EXPECT_NEAR( croty.cos() , cos , eps ) ;
    EXPECT_NEAR( croty.sin() , sin , eps ) ;

    // update(theta)
    croty.update( pi ) ;
    EXPECT_NEAR( croty.theta() , pi , eps ) ;
    EXPECT_NEAR( croty.cos() , 0 , eps ) ;
    EXPECT_NEAR( croty.sin() , 1 , eps ) ;

    // update(cos,sin)
    croty.update( cos , sin ) ;
    EXPECT_NEAR( croty.theta() , pi/2 , eps ) ;
    EXPECT_NEAR( croty.cos() , cos , eps ) ;
    EXPECT_NEAR( croty.sin() , sin , eps ) ;
  }

  {
    const R theta = pi/4 ;
    qclab::QAngle< R >  angle( theta/2 ) ;
    qclab::qgates::CRotationY< T >  croty( 1 , 3 , angle ) ;

    EXPECT_EQ( croty.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( croty.fixed() ) ;          // fixed
    EXPECT_TRUE( croty.controlled() ) ;      // controlled
    EXPECT_EQ( croty.control() , 1 ) ;       // control
    EXPECT_EQ( croty.target() , 3 ) ;        // target
    EXPECT_EQ( croty.controlState() , 1 ) ;  // controlState

    EXPECT_NEAR( croty.theta() , theta , eps ) ;
    EXPECT_NEAR( croty.cos() , std::cos( theta/2 ) , eps ) ;
    EXPECT_NEAR( croty.sin() , std::sin( theta/2 ) , eps ) ;
  }

  {
    const R theta = pi/4 ;
    qclab::qgates::CRotationY< T >  croty( 1 , 3 , theta ) ;

    EXPECT_EQ( croty.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( croty.fixed() ) ;          // fixed
    EXPECT_TRUE( croty.controlled() ) ;      // controlled
    EXPECT_EQ( croty.control() , 1 ) ;       // control
    EXPECT_EQ( croty.target() , 3 ) ;        // target
    EXPECT_EQ( croty.controlState() , 1 ) ;  // controlState

    EXPECT_NEAR( croty.theta() , theta , eps ) ;
    EXPECT_NEAR( croty.cos() , std::cos( theta/2 ) , eps ) ;
    EXPECT_NEAR( croty.sin() , std::sin( theta/2 ) , eps ) ;
  }

  {
    const R theta = pi/4 ;
    const R cos = std::cos( theta ) ;
    const R sin = std::sin( theta ) ;
    qclab::qgates::CRotationY< T >  croty( 1 , 3 , cos , sin ) ;

    EXPECT_EQ( croty.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( croty.fixed() ) ;          // fixed
    EXPECT_TRUE( croty.controlled() ) ;      // controlled
    EXPECT_EQ( croty.control() , 1 ) ;       // control
    EXPECT_EQ( croty.target() , 3 ) ;        // target
    EXPECT_EQ( croty.controlState() , 1 ) ;  // controlState

    EXPECT_NEAR( croty.theta() , 2*theta , eps ) ;
    EXPECT_NEAR( croty.cos() , cos , eps ) ;
    EXPECT_NEAR( croty.sin() , sin , eps ) ;
  }

}


/*
 * complex float
 */
TEST( qclab_qgates_CRotationY , complex_float ) {
  test_qclab_qgates_CRotationY< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_CRotationY , complex_double ) {
  test_qclab_qgates_CRotationY< std::complex< double > >() ;
}

