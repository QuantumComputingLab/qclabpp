#include <gtest/gtest.h>
#include "qgates/RotationZZ.hpp"

template <typename T>
void test_qclab_qgates_RotationZZ() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R eps = std::numeric_limits< R >::epsilon() ;

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

    // update(angle)
    qclab::QAngle< R >  new_angle( 0.5 ) ;
    Rzz.update( new_angle ) ;
    EXPECT_NEAR( Rzz.theta() , 1 , eps ) ;
    EXPECT_NEAR( Rzz.cos() , std::cos( 0.5 ) , eps ) ;
    EXPECT_NEAR( Rzz.sin() , std::sin( 0.5 ) , eps ) ;

    // update(theta)
    Rzz.update( pi/2 ) ;
    EXPECT_NEAR( Rzz.theta() , pi/2 , eps ) ;
    EXPECT_NEAR( Rzz.cos() , std::cos( pi/4 ) , eps ) ;
    EXPECT_NEAR( Rzz.sin() , std::sin( pi/4 ) , eps ) ;

    // matrix
    const T a = std::exp( T(0,-pi/4) ) ;
    const T b = std::exp( T(0, pi/4) ) ;
    EXPECT_NEAR( std::real( Rzz.matrix()(0,0) ) , std::real( a ) , eps ) ;
    EXPECT_NEAR( std::imag( Rzz.matrix()(0,0) ) , std::imag( a ) , eps ) ;
    EXPECT_NEAR( std::real( Rzz.matrix()(1,1) ) , std::real( b ) , eps ) ;
    EXPECT_NEAR( std::imag( Rzz.matrix()(1,1) ) , std::imag( b ) , eps ) ;
    EXPECT_NEAR( std::real( Rzz.matrix()(2,2) ) , std::real( b ) , eps ) ;
    EXPECT_NEAR( std::imag( Rzz.matrix()(2,2) ) , std::imag( b ) , eps ) ;
    EXPECT_NEAR( std::real( Rzz.matrix()(3,3) ) , std::real( a ) , eps ) ;
    EXPECT_NEAR( std::imag( Rzz.matrix()(3,3) ) , std::imag( a ) , eps ) ;
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
    EXPECT_NEAR( Rzz.cos() , std::cos( pi/3 ) , eps ) ;
    EXPECT_NEAR( Rzz.sin() , std::sin( pi/3 ) , eps ) ;
    EXPECT_NEAR( Rzz.theta() , 2*(pi/3) , eps ) ;

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
    qclab::QAngle< R >              angle( pi/4 ) ;
    qclab::qgates::RotationZZ< T >  Rzz( angle ) ;

    EXPECT_EQ( Rzz.nbQubits() , 2 ) ;                    // nbQubits
    EXPECT_FALSE( Rzz.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rzz.controlled() ) ;                   // controlled

    auto qubits = Rzz.qubits() ;
    EXPECT_EQ( qubits[0] , 0 ) ;                         // qubit0
    EXPECT_EQ( qubits[1] , 1 ) ;                         // qubit1

    EXPECT_NEAR( Rzz.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Rzz.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( Rzz.sin() , std::sin( pi/4 ) , eps ) ;  // sin
  }

  {
    qclab::qgates::RotationZZ< T >  Rzz( pi/2 ) ;

    EXPECT_EQ( Rzz.nbQubits() , 2 ) ;                    // nbQubits
    EXPECT_FALSE( Rzz.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rzz.controlled() ) ;                   // controlled

    auto qubits = Rzz.qubits() ;
    EXPECT_EQ( qubits[0] , 0 ) ;                         // qubit0
    EXPECT_EQ( qubits[1] , 1 ) ;                         // qubit1

    EXPECT_NEAR( Rzz.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Rzz.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( Rzz.sin() , std::sin( pi/4 ) , eps ) ;  // sin
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

    EXPECT_NEAR( Rzz.theta() , pi/2 , eps ) ;  // theta
    EXPECT_NEAR( Rzz.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( Rzz.sin() , sin , eps ) ;     // sin
  }


  //
  // constructors with *qubits
  //
  {
    int qbits[2] = { 3 , 5 } ;
    qclab::QAngle< R >              angle( pi/4 ) ;
    qclab::qgates::RotationZZ< T >  Rzz( &qbits[0] , angle ) ;

    EXPECT_EQ( Rzz.nbQubits() , 2 ) ;                    // nbQubits
    EXPECT_FALSE( Rzz.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rzz.controlled() ) ;                   // controlled

    auto qubits = Rzz.qubits() ;
    EXPECT_EQ( qubits[0] , 3 ) ;                         // qubit0
    EXPECT_EQ( qubits[1] , 5 ) ;                         // qubit1

    EXPECT_NEAR( Rzz.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Rzz.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( Rzz.sin() , std::sin( pi/4 ) , eps ) ;  // sin
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

    EXPECT_NEAR( Rzz.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Rzz.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( Rzz.sin() , std::sin( pi/4 ) , eps ) ;  // sin
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

    EXPECT_NEAR( Rzz.theta() , pi/2 , eps ) ;  // theta
    EXPECT_NEAR( Rzz.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( Rzz.sin() , sin , eps ) ;     // sin
  }


  //
  // constructors with qubit0 and qubit1
  //
  {
    int qubit0 = 3 ;
    int qubit1 = 5 ;
    qclab::QAngle< R >              angle( pi/4 ) ;
    qclab::qgates::RotationZZ< T >  Rzz( qubit0 , qubit1 , angle ) ;

    EXPECT_EQ( Rzz.nbQubits() , 2 ) ;                    // nbQubits
    EXPECT_FALSE( Rzz.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rzz.controlled() ) ;                   // controlled

    auto qubits = Rzz.qubits() ;
    EXPECT_EQ( qubits[0] , 3 ) ;                         // qubit0
    EXPECT_EQ( qubits[1] , 5 ) ;                         // qubit1

    EXPECT_NEAR( Rzz.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Rzz.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( Rzz.sin() , std::sin( pi/4 ) , eps ) ;  // sin
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

    EXPECT_NEAR( Rzz.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Rzz.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( Rzz.sin() , std::sin( pi/4 ) , eps ) ;  // sin
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

    EXPECT_NEAR( Rzz.theta() , pi/2 , eps ) ;  // theta
    EXPECT_NEAR( Rzz.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( Rzz.sin() , sin , eps ) ;     // sin
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

