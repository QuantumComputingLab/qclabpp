#include <gtest/gtest.h>
#include "qclab/QCircuit.hpp"
#include "qclab/qgates/Hadamard.hpp"
#include "qclab/qgates/PauliX.hpp"
#include "qclab/qgates/PauliY.hpp"
#include "qclab/qgates/PauliZ.hpp"
#include "qclab/qgates/CNOT.hpp"
#include "qclab/qgates/SWAP.hpp"

template <typename T>
void test_qclab_QCircuit() {

  const T pi = 4 * std::atan(1) ;
  const T eps = std::numeric_limits< T >::epsilon() ;

  using O = std::unique_ptr< qclab::QObject< T > > ;
  using H = qclab::qgates::Hadamard< T > ;
  using X = qclab::qgates::PauliX< T > ;
  using Y = qclab::qgates::PauliY< T > ;
  using Z = qclab::qgates::PauliZ< T > ;
  using CNOT = qclab::qgates::CNOT< T > ;
  using SWAP = qclab::qgates::SWAP< T > ;
  using iter   = typename qclab::QCircuit< T >::iterator ;
  using citer  = typename qclab::QCircuit< T >::const_iterator ;
  using riter  = typename qclab::QCircuit< T >::reverse_iterator ;
  using criter = typename qclab::QCircuit< T >::const_reverse_iterator ;

  {
    qclab::QCircuit< T >  circuit( 1 ) ;

    EXPECT_EQ( circuit.nbQubits() , 1 ) ;     // nbQubits
    EXPECT_FALSE( circuit.fixed() ) ;         // fixed
    EXPECT_FALSE( circuit.controlled() ) ;    // controlled

    // qubit(s)
    EXPECT_EQ( circuit.qubit() , 0 ) ;
    EXPECT_EQ( circuit.qubits().size() , 1 ) ;
    EXPECT_EQ( circuit.qubits()[0] , 0 ) ;
  }

  {
    qclab::QCircuit< T >  circuit( 3 ) ;

    EXPECT_EQ( circuit.nbQubits() , 3 ) ;     // nbQubits
    EXPECT_FALSE( circuit.fixed() ) ;         // fixed
    EXPECT_FALSE( circuit.controlled() ) ;    // controlled

    // qubit(s)
    EXPECT_EQ( circuit.qubit() , 0 ) ;
    EXPECT_EQ( circuit.qubits().size() , 3 ) ;
    EXPECT_EQ( circuit.qubits()[0] , 0 ) ;
    EXPECT_EQ( circuit.qubits()[1] , 1 ) ;
    EXPECT_EQ( circuit.qubits()[2] , 2 ) ;
  }

  {
    qclab::QCircuit< T >  circuit( 3 , 5 ) ;

    EXPECT_EQ( circuit.nbQubits() , 3 ) ;     // nbQubits
    EXPECT_FALSE( circuit.fixed() ) ;         // fixed
    EXPECT_FALSE( circuit.controlled() ) ;    // controlled

    // qubit(s)
    EXPECT_EQ( circuit.qubit() , 5 ) ;
    EXPECT_EQ( circuit.qubits().size() , 3 ) ;
    EXPECT_EQ( circuit.qubits()[0] , 5 ) ;
    EXPECT_EQ( circuit.qubits()[1] , 6 ) ;
    EXPECT_EQ( circuit.qubits()[2] , 7 ) ;

    // offset
    EXPECT_EQ( circuit.offset() , 5 ) ;
    circuit.setOffset( 2 ) ;
    EXPECT_EQ( circuit.offset() , 2 ) ;

    EXPECT_EQ( circuit.qubit() , 2 ) ;
    EXPECT_EQ( circuit.qubits().size() , 3 ) ;
    EXPECT_EQ( circuit.qubits()[0] , 2 ) ;
    EXPECT_EQ( circuit.qubits()[1] , 3 ) ;
    EXPECT_EQ( circuit.qubits()[2] , 4 ) ;
  }

  {
    qclab::QCircuit< T >  circuit( 3 , 0 , 4 ) ;

    EXPECT_EQ( circuit.nbQubits() , 3 ) ;     // nbQubits
    EXPECT_FALSE( circuit.fixed() ) ;         // fixed
    EXPECT_FALSE( circuit.controlled() ) ;    // controlled
    EXPECT_EQ( circuit.nbGates() , 4 ) ;      // nbGates

    // qubit(s)
    EXPECT_EQ( circuit.qubit() , 0 ) ;
    EXPECT_EQ( circuit.qubits().size() , 3 ) ;
    EXPECT_EQ( circuit.qubits()[0] , 0 ) ;
    EXPECT_EQ( circuit.qubits()[1] , 1 ) ;
    EXPECT_EQ( circuit.qubits()[2] , 2 ) ;

    // gates
    EXPECT_TRUE( circuit[0] == nullptr ) ;
    EXPECT_TRUE( circuit[1] == nullptr ) ;
    EXPECT_TRUE( circuit[2] == nullptr ) ;
    EXPECT_TRUE( circuit[3] == nullptr ) ;
  }

  {
    qclab::QCircuit< T , qclab::qgates::QGate1< T > >  circuit( 3 ) ;

    EXPECT_EQ( circuit.nbQubits() , 3 ) ;     // nbQubits
    EXPECT_FALSE( circuit.fixed() ) ;         // fixed
    EXPECT_FALSE( circuit.controlled() ) ;    // controlled

    // qubit(s)
    EXPECT_EQ( circuit.qubit() , 0 ) ;
    EXPECT_EQ( circuit.qubits().size() , 3 ) ;
    EXPECT_EQ( circuit.qubits()[0] , 0 ) ;
    EXPECT_EQ( circuit.qubits()[1] , 1 ) ;
    EXPECT_EQ( circuit.qubits()[2] , 2 ) ;

    // add 1-qubit gates
    circuit.push_back( std::make_unique< X >( 0 ) ) ;
    circuit.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit.push_back( std::make_unique< Z >( 2 ) ) ;
    circuit.push_back( std::make_unique< H >( 1 ) ) ;
    EXPECT_TRUE( *circuit[0] == X( 0 ) ) ;
    EXPECT_TRUE( *circuit[1] == Y( 1 ) ) ;
    EXPECT_TRUE( *circuit[2] == Z( 2 ) ) ;
    EXPECT_TRUE( *circuit[3] == H( 1 ) ) ;
  }

  {
    qclab::QCircuit< T , qclab::qgates::QGate1< T > >  circuit1( 1 ) ;
    qclab::QCircuit< T , qclab::qgates::QGate1< T > >  circuit3( 3 ) ;

    auto I1 = qclab::dense::eye< T >( 2 ) ;
    auto I2 = qclab::dense::eye< T >( 4 ) ;

    // matrix
    circuit1.push_back( std::make_unique< X >( 0 ) ) ;
    EXPECT_TRUE( circuit1.matrix() == X( 0 ).matrix() ) ;
    circuit1.push_back( std::make_unique< Y >( 0 ) ) ;
    EXPECT_TRUE( circuit1.matrix() == Y( 0 ).matrix() * X( 0 ).matrix() ) ;
    circuit1.push_back( std::make_unique< Z >( 0 ) ) ;
    qclab::dense::SquareMatrix< T >  mat1( T(0,-1) ,    0    ,
                                              0    , T(0,-1) ) ;
    EXPECT_TRUE( circuit1.matrix() == mat1 ) ;

    circuit3.push_back( std::make_unique< H >( 0 ) ) ;
    circuit3.push_back( std::make_unique< H >( 1 ) ) ;
    circuit3.push_back( std::make_unique< H >( 2 ) ) ;
    mat1 = H().matrix() ;
    auto mat3 = qclab::dense::kron( mat1 , qclab::dense::kron( mat1 , mat1 ) ) ;
    EXPECT_TRUE( circuit3.matrix() == mat3 ) ;
    circuit3.push_back( std::make_unique< Z >( 0 ) ) ;
    mat3 = qclab::dense::kron( Z().matrix() , I2 ) * mat3 ;
    circuit3.push_back( std::make_unique< Y >( 1 ) ) ;
    mat3 = ( qclab::dense::kron( I1 ,
               qclab::dense::kron( Y().matrix() , I1 ) ) ) * mat3;
    circuit3.push_back( std::make_unique< X >( 2 ) ) ;
    mat3 = qclab::dense::kron( I2 , X().matrix() ) * mat3 ;
    EXPECT_TRUE( circuit3.matrix() == mat3 ) ;
  }

  {
    qclab::QCircuit< T , qclab::qgates::QGate1< T > >  circuit1( 1 ) ;

    auto I1 = qclab::dense::eye< T >( 2 ) ;
    auto I2 = qclab::dense::eye< T >( 4 ) ;

    circuit1.push_back( std::make_unique< H >( 0 ) ) ;
    circuit1.push_back( std::make_unique< Y >( 0 ) ) ;
    circuit1.push_back( std::make_unique< Z >( 0 ) ) ;
    auto mat = I1 ;
    mat = H().matrix() * mat ;
    mat = Y().matrix() * mat ;
    mat = Z().matrix() * mat ;

    // apply (nbQubits = 1)
    auto mat1 = I1 ;
    circuit1.apply( qclab::Side::Left , qclab::Op::NoTrans , 1 , mat1 ) ;
    EXPECT_EQ( mat1 , mat ) ;
    mat1 = I1 ;
    circuit1.apply( qclab::Side::Right , qclab::Op::NoTrans , 1 , mat1 ) ;
    EXPECT_EQ( mat1 , mat ) ;

    // apply (nbQubits = 2)
    auto mat2 = I2 ;
    circuit1.apply( qclab::Side::Left , qclab::Op::NoTrans , 2 , mat2 ) ;
    EXPECT_EQ( mat2 , qclab::dense::kron( mat , I1 ) ) ;
    mat2 = I2 ;
    circuit1.apply( qclab::Side::Right , qclab::Op::NoTrans , 2 , mat2 ) ;
    EXPECT_EQ( mat2 , qclab::dense::kron( mat , I1 ) ) ;

    qclab::dense::conjTransInPlace( mat ) ;

    // apply ConjTrans (nbQubits = 1)
    mat1 = I1 ;
    circuit1.apply( qclab::Side::Left , qclab::Op::ConjTrans , 1 , mat1 ) ;
    EXPECT_EQ( mat1 , mat ) ;
    mat1 = I1 ;
    circuit1.apply( qclab::Side::Right , qclab::Op::ConjTrans , 1 , mat1 ) ;
    EXPECT_EQ( mat1 , mat ) ;

    // apply ConjTrans (nbQubits = 2)
    mat2 = I2 ;
    circuit1.apply( qclab::Side::Left , qclab::Op::ConjTrans , 2 , mat2 ) ;
    EXPECT_EQ( mat2 , qclab::dense::kron( mat , I1 ) ) ;
    mat2 = I2 ;
    circuit1.apply( qclab::Side::Right , qclab::Op::ConjTrans , 2 , mat2 ) ;
    EXPECT_EQ( mat2 , qclab::dense::kron( mat , I1 ) ) ;
  }

  {
    qclab::QCircuit< T >  circuit( 3 ) ;

    circuit.push_back( std::make_unique< X >( 0 ) ) ;
    circuit.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit.push_back( std::make_unique< Z >( 2 ) ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( circuit.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[0];\ny q[1];\nz q[2];\n" ) ;
    std::cout << qasm.str() ;
 }

 {
    qclab::QCircuit< T >  circuit1( 3 ) ;
    circuit1.push_back( std::make_unique< X >( 0 ) ) ;
    circuit1.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit1.push_back( std::make_unique< Z >( 2 ) ) ;

    qclab::QCircuit< T >  circuit2( 3 ) ;
    circuit2.push_back( std::make_unique< Z >( 2 ) ) ;
    circuit2.push_back( std::make_unique< X >( 0 ) ) ;
    circuit2.push_back( std::make_unique< Y >( 1 ) ) ;

    // equals
    EXPECT_TRUE(  circuit1 == circuit2 ) ;
    EXPECT_FALSE( circuit1 != circuit2 ) ;
  }

  //
  // Element access
  //

  {
    qclab::QCircuit< T >  circuit( 3 ) ;

    circuit.push_back( std::make_unique< X >( 0 ) ) ;
    circuit.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit.push_back( std::make_unique< Z >( 2 ) ) ;
    circuit.push_back( std::make_unique< H >( 1 ) ) ;

    // at
    EXPECT_TRUE( *circuit.at( 0 ) == X( 0 ) ) ;
    EXPECT_TRUE( *circuit.at( 0 ) != Y( 1 ) ) ;
    EXPECT_TRUE( *circuit.at( 0 ) != Z( 2 ) ) ;
    EXPECT_TRUE( *circuit.at( 1 ) == Y( 1 ) ) ;
    EXPECT_TRUE( *circuit.at( 2 ) == Z( 2 ) ) ;
    EXPECT_TRUE( *circuit.at( 3 ) == H( 1 ) ) ;
    {
      O& x = circuit.at( 0 ) ; EXPECT_TRUE( *x == X( 0 ) ) ;
      O& y = circuit.at( 1 ) ; EXPECT_TRUE( *y == Y( 1 ) ) ;
      O& z = circuit.at( 2 ) ; EXPECT_TRUE( *z == Z( 2 ) ) ;
      O& h = circuit.at( 3 ) ; EXPECT_TRUE( *h == H( 1 ) ) ;
    }

    // operator[]
    EXPECT_TRUE( *circuit[0] == X( 0 ) ) ;
    EXPECT_TRUE( *circuit[1] == Y( 1 ) ) ;
    EXPECT_TRUE( *circuit[2] == Z( 2 ) ) ;
    EXPECT_TRUE( *circuit[3] == H( 1 ) ) ;
    {
      O& x = circuit[0] ; EXPECT_TRUE( *x == X( 0 ) ) ;
      O& y = circuit[1] ; EXPECT_TRUE( *y == Y( 1 ) ) ;
      O& z = circuit[2] ; EXPECT_TRUE( *z == Z( 2 ) ) ;
      O& h = circuit[3] ; EXPECT_TRUE( *h == H( 1 ) ) ;
    }

    // front
    EXPECT_TRUE( *circuit.front() == X(0) ) ;
    O& x = circuit.front() ; EXPECT_TRUE( *x == X() ) ;

    // back
    EXPECT_TRUE( *circuit.back() == H(1) ) ;
    O& h = circuit.back() ; EXPECT_TRUE( *h == H() ) ;
  }

  //
  // Iterators
  //

  {
    qclab::QCircuit< T >  circuit( 3 ) ;

    circuit.push_back( std::make_unique< X >( 0 ) ) ;
    circuit.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit.push_back( std::make_unique< Z >( 2 ) ) ;
    circuit.push_back( std::make_unique< H >( 1 ) ) ;

    // begin
    {
      EXPECT_TRUE( **circuit.begin() == X() ) ;
      iter it1 = circuit.begin() ;
      EXPECT_TRUE( **it1 == X() ) ;
      citer it2 = circuit.begin() ;
      EXPECT_TRUE( **it2 == X() ) ;
    }

    // end
    {
      EXPECT_TRUE( circuit.end() == circuit.begin() + 4 ) ;
      iter it1 = circuit.end() ;
      EXPECT_TRUE( it1 == circuit.begin() + 4 ) ;
      citer it2 = circuit.end() ;
      EXPECT_TRUE( it2 == circuit.begin() + 4 ) ;
    }

    // rbegin
    {
      EXPECT_TRUE( **circuit.rbegin() == H() ) ;
      riter it1 = circuit.rbegin() ;
      EXPECT_TRUE( **it1 == H() ) ;
      criter it2 = circuit.rbegin() ;
      EXPECT_TRUE( **it2 == H() ) ;
    }

    // rend
    {
      EXPECT_TRUE( circuit.rend() == circuit.rbegin() + 4 ) ;
      riter it1 = circuit.rend() ;
      EXPECT_TRUE( it1 == circuit.rbegin() + 4 ) ;
      criter it2 = circuit.rend() ;
      EXPECT_TRUE( it2 == circuit.rbegin() + 4 ) ;
    }
  }

  //
  // Capacity
  //

  {
    qclab::QCircuit< T >  circuit( 3 ) ;

    // empty & nbGates
    EXPECT_TRUE( circuit.empty() ) ;
    EXPECT_EQ( circuit.nbGates() , 0 ) ;
    circuit.push_back( std::make_unique< X >( 0 ) ) ;
    EXPECT_FALSE( circuit.empty() ) ;
    EXPECT_EQ( circuit.nbGates() , 1 ) ;

    // reserve & capacity
    EXPECT_TRUE( circuit.capacity() != 5 ) ;
    circuit.reserve( 5 ) ;
    EXPECT_EQ( circuit.capacity() , 5 ) ;

    // shrink_to_fit
    circuit.shrink_to_fit() ;
    EXPECT_EQ( circuit.capacity() , 1 ) ;
  }

  //
  // Modifiers
  //

  {
    qclab::QCircuit< T >  circuit( 3 ) ;

    circuit.push_back( std::make_unique< X >( 0 ) ) ;
    circuit.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit.push_back( std::make_unique< Z >( 2 ) ) ;
    circuit.push_back( std::make_unique< H >( 1 ) ) ;

    // clear
    EXPECT_FALSE( circuit.empty() ) ;
    EXPECT_EQ( circuit.nbGates() , 4 ) ;
    circuit.clear() ;
    EXPECT_TRUE( circuit.empty() ) ;
    EXPECT_EQ( circuit.nbGates() , 0 ) ;
  }

  {
    qclab::QCircuit< T >  circuit( 3 ) ;

    // insert 1 gate
    auto it1 = circuit.insert( circuit.begin() , std::make_unique< X >( 0 ) ) ;
    EXPECT_EQ( circuit.nbGates() , 1 ) ;
    EXPECT_TRUE( *circuit[0] == X( 0 ) ) ;
    EXPECT_TRUE( it1 == circuit.begin() ) ;
    EXPECT_TRUE( **it1 == X( 0 ) ) ;
    auto it2 = circuit.insert( circuit.begin() , std::make_unique< Y >( 1 ) ) ;
    EXPECT_EQ( circuit.nbGates() , 2 ) ;
    EXPECT_TRUE( *circuit[0] == Y( 1 ) ) ;
    EXPECT_TRUE( *circuit[1] == X( 0 ) ) ;
    EXPECT_TRUE( it2 == circuit.begin() ) ;
    EXPECT_TRUE( **it2 == Y( 1 ) ) ;
    auto it3 = circuit.insert( circuit.end() , std::make_unique< Z >( 2 ) ) ;
    EXPECT_EQ( circuit.nbGates() , 3 ) ;
    EXPECT_TRUE( *circuit[0] == Y( 1 ) ) ;
    EXPECT_TRUE( *circuit[1] == X( 0 ) ) ;
    EXPECT_TRUE( *circuit[2] == Z( 2 ) ) ;
    EXPECT_TRUE( it3 == circuit.end() - 1 ) ;
    EXPECT_TRUE( **it3 == Z( 2 ) ) ;
    auto it4 = circuit.insert( circuit.begin() + 1 , std::make_unique< H >() ) ;
    EXPECT_EQ( circuit.nbGates() , 4 ) ;
    EXPECT_TRUE( *circuit[0] == Y( 1 ) ) ;
    EXPECT_TRUE( *circuit[1] == H( 0 ) ) ;
    EXPECT_TRUE( *circuit[2] == X( 0 ) ) ;
    EXPECT_TRUE( *circuit[3] == Z( 2 ) ) ;
    EXPECT_TRUE( it4 == circuit.begin() + 1 ) ;
    EXPECT_TRUE( **it4 == H( 0 ) ) ;

    EXPECT_DEBUG_DEATH( circuit.insert( circuit.begin() ,
                                        std::make_unique< H >( 3 ) ) , "" ) ;
    EXPECT_DEBUG_DEATH( circuit.insert( circuit.begin() ,
                                        std::make_unique< H >( 4 ) ) , "" ) ;

    circuit.insert( circuit.end() , std::make_unique< CNOT >( 0 , 1 ) ) ;
    circuit.insert( circuit.end() , std::make_unique< CNOT >( 1 , 2 ) ) ;
    EXPECT_DEBUG_DEATH( circuit.insert( circuit.end() ,
                                        std::make_unique< CNOT >( 2 , 3 ) ) ,
                        "" ) ;
    EXPECT_DEBUG_DEATH( circuit.insert( circuit.end() ,
                                        std::make_unique< CNOT >( 3 , 4 ) ) ,
                        "" ) ;
  }

  {
    qclab::QCircuit< T >  circuit1( 3 ) ;

    circuit1.push_back( std::make_unique< X >( 0 ) ) ;
    circuit1.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit1.push_back( std::make_unique< Z >( 2 ) ) ;

    // insert all
    qclab::QCircuit< T >  circuit2( 3 ) ;
    auto begin1 = circuit1.begin() ;
    auto end1 = circuit1.end() ;
    iter it2 = circuit2.insert( circuit2.begin() , begin1 , end1 ) ;
    EXPECT_EQ( circuit2.nbGates() , 3 ) ;
    EXPECT_TRUE( *circuit2[0] == X( 0 ) ) ;
    EXPECT_TRUE( *circuit2[1] == Y( 1 ) ) ;
    EXPECT_TRUE( *circuit2[2] == Z( 2 ) ) ;
    EXPECT_EQ( circuit1[0] , nullptr ) ;
    EXPECT_EQ( circuit1[1] , nullptr ) ;
    EXPECT_EQ( circuit1[2] , nullptr ) ;
    EXPECT_TRUE( it2 == circuit2.begin() ) ;
    EXPECT_TRUE( **it2 == X( 0 ) ) ;
  }

  {
    qclab::QCircuit< T >  circuit1( 3 ) ;

    circuit1.push_back( std::make_unique< X >( 0 ) ) ;
    circuit1.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit1.push_back( std::make_unique< Z >( 2 ) ) ;

    // insert all before
    qclab::QCircuit< T >  circuit2( 3 ) ;
    circuit2.push_back( std::make_unique< H >( 1 ) ) ;
    auto begin1 = circuit1.begin() ;
    auto end1 = circuit1.end() ;
    iter it2 = circuit2.insert( circuit2.begin() , begin1 , end1 ) ;
    EXPECT_EQ( circuit2.nbGates() , 4 ) ;
    EXPECT_TRUE( *circuit2[0] == X( 0 ) ) ;
    EXPECT_TRUE( *circuit2[1] == Y( 1 ) ) ;
    EXPECT_TRUE( *circuit2[2] == Z( 2 ) ) ;
    EXPECT_EQ( circuit1[0] , nullptr ) ;
    EXPECT_EQ( circuit1[1] , nullptr ) ;
    EXPECT_EQ( circuit1[2] , nullptr ) ;
    EXPECT_TRUE( it2 == circuit2.begin() ) ;
    EXPECT_TRUE( **it2 == X( 0 ) ) ;
  }

  {
    qclab::QCircuit< T >  circuit1( 3 ) ;

    circuit1.push_back( std::make_unique< X >( 0 ) ) ;
    circuit1.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit1.push_back( std::make_unique< Z >( 2 ) ) ;

    // insert all after
    qclab::QCircuit< T >  circuit2( 3 ) ;
    circuit2.push_back( std::make_unique< H >( 1 ) ) ;
    auto begin1 = circuit1.begin() ;
    auto end1 = circuit1.end() ;
    iter it2 = circuit2.insert( circuit2.end() , begin1 , end1 ) ;
    EXPECT_EQ( circuit2.nbGates() , 4 ) ;
    EXPECT_TRUE( *circuit2[1] == X( 0 ) ) ;
    EXPECT_TRUE( *circuit2[2] == Y( 1 ) ) ;
    EXPECT_TRUE( *circuit2[3] == Z( 2 ) ) ;
    EXPECT_EQ( circuit1[0] , nullptr ) ;
    EXPECT_EQ( circuit1[1] , nullptr ) ;
    EXPECT_EQ( circuit1[2] , nullptr ) ;
    EXPECT_TRUE( it2 == circuit2.begin() + 1 ) ;
    EXPECT_TRUE( **it2 == X( 0 ) ) ;
  }

  {
    qclab::QCircuit< T >  circuit1( 3 ) ;

    circuit1.push_back( std::make_unique< X >( 0 ) ) ;
    circuit1.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit1.push_back( std::make_unique< Z >( 2 ) ) ;

    // insert nothing before
    qclab::QCircuit< T >  circuit2( 3 ) ;
    circuit2.push_back( std::make_unique< H >( 1 ) ) ;
    auto begin1 = circuit1.begin() ;
    iter it2 = circuit2.insert( circuit2.begin() , begin1 , begin1 ) ;
    EXPECT_EQ( circuit2.nbGates() , 1 ) ;
    EXPECT_TRUE( it2 == circuit2.begin() ) ;
  }

  {
    qclab::QCircuit< T >  circuit1( 3 ) ;

    circuit1.push_back( std::make_unique< X >( 0 ) ) ;
    circuit1.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit1.push_back( std::make_unique< Z >( 2 ) ) ;

    // insert nothing after
    qclab::QCircuit< T >  circuit2( 3 ) ;
    circuit2.push_back( std::make_unique< H >( 1 ) ) ;
    auto begin1 = circuit1.begin() ;
    iter it2 = circuit2.insert( circuit2.end() , begin1 , begin1 ) ;
    EXPECT_EQ( circuit2.nbGates() , 1 ) ;
    EXPECT_TRUE( it2 == circuit2.end() ) ;
  }

  {
    qclab::QCircuit< T >  circuit( 3 ) ;

    circuit.push_back( std::make_unique< X >( 0 ) ) ;
    circuit.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit.push_back( std::make_unique< Z >( 2 ) ) ;

    // erase
    auto it1 = circuit.erase( circuit.begin() + 1 ) ;
    EXPECT_EQ( circuit.nbGates() , 2 ) ;
    EXPECT_TRUE( *circuit[0] == X( 0 ) ) ;
    EXPECT_TRUE( *circuit[1] == Z( 2 ) ) ;
    EXPECT_TRUE( it1 == circuit.begin() + 1 ) ;

    auto it2 = circuit.erase( circuit.begin() + 1 ) ;
    EXPECT_EQ( circuit.nbGates() , 1 ) ;
    EXPECT_TRUE( *circuit[0] == X( 0 ) ) ;
    EXPECT_TRUE( it2 == circuit.end() ) ;

    auto it3 = circuit.erase( circuit.begin() ) ;
    EXPECT_EQ( circuit.nbGates() , 0 ) ;
    EXPECT_TRUE( circuit.empty() ) ;
    EXPECT_TRUE( it3 == circuit.end() ) ;
  }

  {
    qclab::QCircuit< T >  circuit( 3 ) ;

    circuit.push_back( std::make_unique< X >( 0 ) ) ;
    circuit.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit.push_back( std::make_unique< Z >( 2 ) ) ;

    // erase all
    auto it = circuit.erase( circuit.begin() , circuit.end() ) ;
    EXPECT_TRUE( circuit.empty() ) ;
    EXPECT_TRUE( it == circuit.end() ) ;
  }

  {
    qclab::QCircuit< T >  circuit( 3 ) ;

    circuit.push_back( std::make_unique< X >( 0 ) ) ;
    circuit.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit.push_back( std::make_unique< Z >( 2 ) ) ;

    // erase all before
    auto it = circuit.erase( circuit.begin() , circuit.begin() + 2 ) ;
    EXPECT_EQ( circuit.nbGates() , 1 ) ;
    EXPECT_TRUE( *circuit[0] == Z( 2 ) ) ;
    EXPECT_TRUE( it == circuit.begin() ) ;
  }

  {
    qclab::QCircuit< T >  circuit( 3 ) ;

    circuit.push_back( std::make_unique< X >( 0 ) ) ;
    circuit.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit.push_back( std::make_unique< Z >( 2 ) ) ;

    // erase all after
    auto it = circuit.erase( circuit.begin() + 1 , circuit.end() ) ;
    EXPECT_EQ( circuit.nbGates() , 1 ) ;
    EXPECT_TRUE( *circuit[0] == X( 0 ) ) ;
    EXPECT_TRUE( it == circuit.end() ) ;
  }

  {
    qclab::QCircuit< T >  circuit( 3 ) ;

    circuit.push_back( std::make_unique< X >( 0 ) ) ;
    circuit.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit.push_back( std::make_unique< Z >( 2 ) ) ;

    // erase nothing
    auto it = circuit.erase( circuit.begin() + 1 , circuit.begin() + 1 ) ;
    EXPECT_EQ( circuit.nbGates() , 3 ) ;
    EXPECT_TRUE( *circuit[0] == X( 0 ) ) ;
    EXPECT_TRUE( *circuit[1] == Y( 1 ) ) ;
    EXPECT_TRUE( *circuit[2] == Z( 2 ) ) ;
    EXPECT_TRUE( it == circuit.begin() + 1 ) ;
  }

  {
    qclab::QCircuit< T >  circuit( 3 ) ;

    // push_back
    circuit.push_back( std::make_unique< X >( 0 ) ) ;
    EXPECT_EQ( circuit.nbGates() , 1 ) ;
    EXPECT_TRUE( *circuit[0] == X( 0 ) ) ;
    circuit.push_back( std::make_unique< Y >( 1 ) ) ;
    EXPECT_EQ( circuit.nbGates() , 2 ) ;
    EXPECT_TRUE( *circuit[0] == X( 0 ) ) ;
    EXPECT_TRUE( *circuit[1] == Y( 1 ) ) ;
    circuit.push_back( std::make_unique< Z >( 2 ) ) ;
    EXPECT_EQ( circuit.nbGates() , 3 ) ;
    EXPECT_TRUE( *circuit[0] == X( 0 ) ) ;
    EXPECT_TRUE( *circuit[1] == Y( 1 ) ) ;
    EXPECT_TRUE( *circuit[2] == Z( 2 ) ) ;
  }

  {
    qclab::QCircuit< T >  circuit( 3 ) ;

    circuit.push_back( std::make_unique< X >( 0 ) ) ;
    circuit.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit.push_back( std::make_unique< Z >( 2 ) ) ;

    // pop_back
    circuit.pop_back() ;
    EXPECT_EQ( circuit.nbGates() , 2 ) ;
    EXPECT_TRUE( *circuit[0] == X( 0 ) ) ;
    EXPECT_TRUE( *circuit[1] == Y( 1 ) ) ;
    circuit.pop_back() ;
    EXPECT_EQ( circuit.nbGates() , 1 ) ;
    EXPECT_TRUE( *circuit[0] == X( 0 ) ) ;
    circuit.pop_back() ;
    EXPECT_TRUE( circuit.empty() ) ;
  }

  {
    qclab::QCircuit< T >  circuit( 3 ) ;

    circuit.push_back( std::make_unique< X >( 0 ) ) ;
    circuit.push_back( std::make_unique< Y >( 1 ) ) ;
    circuit.push_back( std::make_unique< Z >( 2 ) ) ;

    // resize
    circuit.resize( 5 ) ;
    EXPECT_EQ( circuit.nbGates() , 5 ) ;
    EXPECT_TRUE( *circuit[0] == X( 0 ) ) ;
    EXPECT_TRUE( *circuit[1] == Y( 1 ) ) ;
    EXPECT_TRUE( *circuit[2] == Z( 2 ) ) ;
    EXPECT_TRUE(  circuit[3] == nullptr ) ;
    EXPECT_TRUE(  circuit[4] == nullptr ) ;
    circuit.resize( 3 ) ;
    EXPECT_EQ( circuit.nbGates() , 3 ) ;
    EXPECT_TRUE( *circuit[0] == X( 0 ) ) ;
    EXPECT_TRUE( *circuit[1] == Y( 1 ) ) ;
    EXPECT_TRUE( *circuit[2] == Z( 2 ) ) ;
    circuit.resize( 2 ) ;
    EXPECT_EQ( circuit.nbGates() , 2 ) ;
    EXPECT_TRUE( *circuit[0] == X( 0 ) ) ;
    EXPECT_TRUE( *circuit[1] == Y( 1 ) ) ;
    circuit.resize( 0 ) ;
    EXPECT_TRUE( circuit.empty() ) ;
  }

}


/*
 * complex float
 */
TEST( qclab_QCircuit , complex_float ) {
  test_qclab_QCircuit< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_QCircuit , complex_double ) {
  test_qclab_QCircuit< std::complex< double > >() ;
}

