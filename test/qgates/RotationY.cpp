#include <gtest/gtest.h>
#include "qgates/RotationY.hpp"

template <typename T>
void test_qclab_qgates_RotationY() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  {
    qclab::qgates::RotationY< T >  Ry ;

    EXPECT_EQ( Ry.nbQubits() , 1 ) ;   // nbQubits
    EXPECT_FALSE( Ry.fixed() ) ;       // fixed
    EXPECT_FALSE( Ry.controlled() ) ;  // controlled
    EXPECT_EQ( Ry.cos() , 1.0 ) ;      // cos
    EXPECT_EQ( Ry.sin() , 0.0 ) ;      // sin
    EXPECT_EQ( Ry.theta() , 0.0 ) ;    // theta

    // matrix
    auto eye = qclab::dense::eye< T >( 2 ) ;
    EXPECT_TRUE( Ry.matrix() == eye ) ;

    // qubit
    EXPECT_EQ( Ry.qubit() , 0 ) ;
    Ry.setQubit( 2 ) ;
    EXPECT_EQ( Ry.qubit() , 2 ) ;

    // qubits
    auto qubits = Ry.qubits() ;
    EXPECT_EQ( qubits.size() , 1 ) ;
    EXPECT_EQ( qubits[0] , 2 ) ;
    int qnew = 3 ;
    Ry.setQubits( &qnew ) ;
    EXPECT_EQ( Ry.qubit() , 3 ) ;

    // fixed
    Ry.makeFixed() ;
    EXPECT_TRUE( Ry.fixed() ) ;
    Ry.makeVariable() ;
    EXPECT_FALSE( Ry.fixed() ) ;

    // update(angle)
    qclab::QAngle< R >  new_angle( 0.5 ) ;
    Ry.update( new_angle ) ;
    EXPECT_NEAR( Ry.theta() , 1 , eps ) ;
    EXPECT_NEAR( Ry.cos() , std::cos( 0.5 ) , eps ) ;
    EXPECT_NEAR( Ry.sin() , std::sin( 0.5 ) , eps ) ;

    // update(theta)
    Ry.update( pi/2 ) ;
    EXPECT_NEAR( Ry.theta() , pi/2 , eps ) ;
    EXPECT_NEAR( Ry.cos() , std::cos( pi/4 ) , eps ) ;
    EXPECT_NEAR( Ry.sin() , std::sin( pi/4 ) , eps ) ;

    // matrix
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    EXPECT_EQ( Ry.matrix()(0,0) , T( cos) ) ;
    EXPECT_EQ( Ry.matrix()(1,0) , T( sin) ) ;
    EXPECT_EQ( Ry.matrix()(0,1) , T(-sin) ) ;
    EXPECT_EQ( Ry.matrix()(1,1) , T( cos) ) ;

    // print
    Ry.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( Ry.toQASM( qasm ) , 0 ) ;
    std::string qasm_check = "ry(" + qclab::qasm( Ry.theta() ) + ") q[3];\n" ;
    EXPECT_EQ( qasm.str() , qasm_check ) ;
    std::cout << qasm.str() ;

    // update(cos,sin)
    Ry.update( std::cos( pi/3 ) , std::sin( pi/3 ) ) ;
    EXPECT_EQ( Ry.cos() , std::cos( pi/3 ) ) ;
    EXPECT_EQ( Ry.sin() , std::sin( pi/3 ) ) ;
    EXPECT_NEAR( Ry.theta() , 2*(pi/3) , eps ) ;

    // operators == and !=
    qclab::qgates::RotationY< T >  Ry2( std::cos( pi/3 ) , std::sin( pi/3 ) ) ;
    EXPECT_TRUE( Ry == Ry2 ) ;
    EXPECT_FALSE( Ry != Ry2 ) ;
    Ry2.update( 1 ) ;
    EXPECT_TRUE( Ry != Ry2 ) ;
    EXPECT_FALSE( Ry == Ry2 ) ;
  }

  {
    qclab::qgates::RotationY< T >  Ry( pi/2 ) ;

    EXPECT_EQ( Ry.nbQubits() , 1 ) ;                    // nbQubits
    EXPECT_FALSE( Ry.fixed() ) ;                        // fixed
    EXPECT_FALSE( Ry.controlled() ) ;                   // controlled
    EXPECT_EQ( Ry.qubit() , 0 ) ;                       // qubit

    EXPECT_NEAR( Ry.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Ry.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( Ry.sin() , std::sin( pi/4 ) , eps ) ;  // sin
  }

  {
    qclab::qgates::RotationY< T >  Ry( 5 , pi/2 ) ;

    EXPECT_EQ( Ry.nbQubits() , 1 ) ;                    // nbQubits
    EXPECT_FALSE( Ry.fixed() ) ;                        // fixed
    EXPECT_FALSE( Ry.controlled() ) ;                   // controlled
    EXPECT_EQ( Ry.qubit() , 5 ) ;                       // qubit

    EXPECT_NEAR( Ry.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Ry.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( Ry.sin() , std::sin( pi/4 ) , eps ) ;  // sin
  }

  {
    qclab::qgates::RotationY< T >  Ry( 5 , pi/2 , true ) ;

    EXPECT_EQ( Ry.nbQubits() , 1 ) ;                    // nbQubits
    EXPECT_TRUE( Ry.fixed() ) ;                         // fixed
    EXPECT_FALSE( Ry.controlled() ) ;                   // controlled
    EXPECT_EQ( Ry.qubit() , 5 ) ;                       // qubit

    EXPECT_NEAR( Ry.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Ry.cos() , std::cos( pi/4 ) , eps ) ;  // cos
    EXPECT_NEAR( Ry.sin() , std::sin( pi/4 ) , eps ) ;  // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationY< T >  Ry( cos , sin ) ;

    EXPECT_EQ( Ry.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_FALSE( Ry.fixed() ) ;              // fixed
    EXPECT_FALSE( Ry.controlled() ) ;         // controlled
    EXPECT_EQ( Ry.qubit() , 0 ) ;             // qubit

    EXPECT_NEAR( Ry.theta() , pi/2 , eps ) ;  // theta
    EXPECT_NEAR( Ry.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( Ry.sin() , sin , eps ) ;     // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationY< T >  Ry( 5 , cos , sin ) ;

    EXPECT_EQ( Ry.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_FALSE( Ry.fixed() ) ;              // fixed
    EXPECT_FALSE( Ry.controlled() ) ;         // controlled
    EXPECT_EQ( Ry.qubit() , 5 ) ;             // qubit

    EXPECT_NEAR( Ry.theta() , pi/2 , eps ) ;  // theta
    EXPECT_NEAR( Ry.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( Ry.sin() , sin , eps ) ;     // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationY< T >  Ry( 5 , cos , sin , true ) ;

    EXPECT_EQ( Ry.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_TRUE( Ry.fixed() ) ;               // fixed
    EXPECT_FALSE( Ry.controlled() ) ;         // controlled
    EXPECT_EQ( Ry.qubit() , 5 ) ;             // qubit

    EXPECT_NEAR( Ry.theta() , pi/2 , eps ) ;  // theta
    EXPECT_NEAR( Ry.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( Ry.sin() , sin , eps ) ;     // sin
  }

}


/*
 * float
 */
TEST( qclab_qgates_RotationY , float ) {
  test_qclab_qgates_RotationY< float >() ;
}

/*
 * double
 */
TEST( qclab_qgates_RotationY , double ) {
  test_qclab_qgates_RotationY< double >() ;
}

/*
 * complex float
 */
TEST( qclab_qgates_RotationY , complex_float ) {
  test_qclab_qgates_RotationY< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_RotationY , complex_double ) {
  test_qclab_qgates_RotationY< std::complex< double > >() ;
}

