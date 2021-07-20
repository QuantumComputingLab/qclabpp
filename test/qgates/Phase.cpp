#include <gtest/gtest.h>
#include "qclab/qgates/Phase.hpp"

template <typename T>
void test_qclab_qgates_Phase() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R eps = std::numeric_limits< R >::epsilon() ;

  {
    qclab::qgates::Phase< T >  Ph ;

    EXPECT_EQ( Ph.nbQubits() , 1 ) ;   // nbQubits
    EXPECT_FALSE( Ph.fixed() ) ;       // fixed
    EXPECT_FALSE( Ph.controlled() ) ;  // controlled
    EXPECT_EQ( Ph.cos() , 1.0 ) ;      // cos
    EXPECT_EQ( Ph.sin() , 0.0 ) ;      // sin
    EXPECT_EQ( Ph.theta() , 0.0 ) ;    // theta

    // matrix
    auto eye = qclab::dense::eye< T >( 2 ) ;
    EXPECT_TRUE( Ph.matrix() == eye ) ;

    // qubit
    EXPECT_EQ( Ph.qubit() , 0 ) ;
    Ph.setQubit( 2 ) ;
    EXPECT_EQ( Ph.qubit() , 2 ) ;

    // qubits
    auto qubits = Ph.qubits() ;
    EXPECT_EQ( qubits.size() , 1 ) ;
    EXPECT_EQ( qubits[0] , 2 ) ;
    int qnew = 3 ;
    Ph.setQubits( &qnew ) ;
    EXPECT_EQ( Ph.qubit() , 3 ) ;

    // fixed
    Ph.makeFixed() ;
    EXPECT_TRUE( Ph.fixed() ) ;
    Ph.makeVariable() ;
    EXPECT_FALSE( Ph.fixed() ) ;

    // update(angle)
    qclab::QAngle< R >  new_angle( 0.5 ) ;
    Ph.update( new_angle ) ;
    EXPECT_NEAR( Ph.theta() , 0.5 , eps ) ;
    EXPECT_NEAR( Ph.cos() , std::cos( 0.5 ) , eps ) ;
    EXPECT_NEAR( Ph.sin() , std::sin( 0.5 ) , eps ) ;

    // update(theta)
    Ph.update( pi/4 ) ;
    EXPECT_NEAR( Ph.theta() , pi/4 , eps ) ;
    EXPECT_NEAR( Ph.cos() , std::cos( pi/4 ) , eps ) ;
    EXPECT_NEAR( Ph.sin() , std::sin( pi/4 ) , eps ) ;

    // matrix
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    EXPECT_EQ( Ph.matrix()(0,0) , T(1) ) ;
    EXPECT_EQ( Ph.matrix()(1,0) , T(0) ) ;
    EXPECT_EQ( Ph.matrix()(0,1) , T(0) ) ;
    EXPECT_EQ( Ph.matrix()(1,1) , T(cos,sin) ) ;

    // print
    Ph.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( Ph.toQASM( qasm ) , 0 ) ;
    std::string qasm_check = "rz(" + qclab::qasm( Ph.theta() ) + ") q[3];\n" ;
    EXPECT_EQ( qasm.str() , qasm_check ) ;
    std::cout << qasm.str() ;

    // update(cos,sin)
    Ph.update( std::cos( pi/3 ) , std::sin( pi/3 ) ) ;
    EXPECT_NEAR( Ph.cos() , std::cos( pi/3 ) , eps ) ;
    EXPECT_NEAR( Ph.sin() , std::sin( pi/3 ) , eps ) ;
    EXPECT_NEAR( Ph.theta() , pi/3 , eps ) ;

    // operators == and !=
    qclab::qgates::Phase< T >  Ph2( std::cos( pi/3 ) , std::sin( pi/3 ) ) ;
    EXPECT_TRUE( Ph == Ph2 ) ;
    EXPECT_FALSE( Ph != Ph2 ) ;
    Ph2.update( 1 ) ;
    EXPECT_TRUE( Ph != Ph2 ) ;
    EXPECT_FALSE( Ph == Ph2 ) ;
  }

  {
    qclab::qgates::Phase< T >  Ph( pi/2 ) ;

    EXPECT_EQ( Ph.nbQubits() , 1 ) ;                    // nbQubits
    EXPECT_FALSE( Ph.fixed() ) ;                        // fixed
    EXPECT_FALSE( Ph.controlled() ) ;                   // controlled
    EXPECT_EQ( Ph.qubit() , 0 ) ;                       // qubit

    EXPECT_NEAR( Ph.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Ph.cos() , std::cos( pi/2 ) , eps ) ;  // cos
    EXPECT_NEAR( Ph.sin() , std::sin( pi/2 ) , eps ) ;  // sin
  }

  {
    qclab::qgates::Phase< T >  Ph( 5 , pi/2 ) ;

    EXPECT_EQ( Ph.nbQubits() , 1 ) ;                    // nbQubits
    EXPECT_FALSE( Ph.fixed() ) ;                        // fixed
    EXPECT_FALSE( Ph.controlled() ) ;                   // controlled
    EXPECT_EQ( Ph.qubit() , 5 ) ;                       // qubit

    EXPECT_NEAR( Ph.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Ph.cos() , std::cos( pi/2 ) , eps ) ;  // cos
    EXPECT_NEAR( Ph.sin() , std::sin( pi/2 ) , eps ) ;  // sin
  }

  {
    qclab::qgates::Phase< T >  Ph( 5 , pi/2 , true ) ;

    EXPECT_EQ( Ph.nbQubits() , 1 ) ;                    // nbQubits
    EXPECT_TRUE( Ph.fixed() ) ;                         // fixed
    EXPECT_FALSE( Ph.controlled() ) ;                   // controlled
    EXPECT_EQ( Ph.qubit() , 5 ) ;                       // qubit

    EXPECT_NEAR( Ph.theta() , pi/2 , eps ) ;            // theta
    EXPECT_NEAR( Ph.cos() , std::cos( pi/2 ) , eps ) ;  // cos
    EXPECT_NEAR( Ph.sin() , std::sin( pi/2 ) , eps ) ;  // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::Phase< T >  Ph( cos , sin ) ;

    EXPECT_EQ( Ph.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_FALSE( Ph.fixed() ) ;              // fixed
    EXPECT_FALSE( Ph.controlled() ) ;         // controlled
    EXPECT_EQ( Ph.qubit() , 0 ) ;             // qubit

    EXPECT_NEAR( Ph.theta() , pi/4 , eps ) ;  // theta
    EXPECT_NEAR( Ph.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( Ph.sin() , sin , eps ) ;     // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::Phase< T >  Ph( 5 , cos , sin ) ;

    EXPECT_EQ( Ph.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_FALSE( Ph.fixed() ) ;              // fixed
    EXPECT_FALSE( Ph.controlled() ) ;         // controlled
    EXPECT_EQ( Ph.qubit() , 5 ) ;             // qubit

    EXPECT_NEAR( Ph.theta() , pi/4 , eps ) ;  // theta
    EXPECT_NEAR( Ph.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( Ph.sin() , sin , eps ) ;     // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::Phase< T >  Ph( 5 , cos , sin , true ) ;

    EXPECT_EQ( Ph.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_TRUE( Ph.fixed() ) ;               // fixed
    EXPECT_FALSE( Ph.controlled() ) ;         // controlled
    EXPECT_EQ( Ph.qubit() , 5 ) ;             // qubit

    EXPECT_NEAR( Ph.theta() , pi/4 , eps ) ;  // theta
    EXPECT_NEAR( Ph.cos() , cos , eps ) ;     // cos
    EXPECT_NEAR( Ph.sin() , sin , eps ) ;     // sin
  }

}


/*
 * complex float
 */
TEST( qclab_qgates_Phase , complex_float ) {
  test_qclab_qgates_Phase< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_Phase , complex_double ) {
  test_qclab_qgates_Phase< std::complex< double > >() ;
}

