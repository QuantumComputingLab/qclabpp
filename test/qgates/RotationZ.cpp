#include <gtest/gtest.h>
#include "qgates/RotationZ.hpp"

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

    // update(angle)
    qclab::QAngle< R >  new_angle( 0.5 ) ;
    Rz.update( new_angle ) ;
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

