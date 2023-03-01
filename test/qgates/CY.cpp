#include <gtest/gtest.h>
#include "qclab/qgates/CY.hpp"

template <typename T>
void test_qclab_qgates_CY() {

  //
  // CYs |1> controlled
  //
  {
    qclab::qgates::CY< T >  cy ;

    EXPECT_EQ( cy.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cy.fixed() ) ;           // fixed
    EXPECT_TRUE( cy.controlled() ) ;      // controlled
    EXPECT_EQ( cy.control() , 0 ) ;       // control
    EXPECT_EQ( cy.target() , 1 ) ;        // target
    EXPECT_EQ( cy.controlState() , 1 ) ;  // controlState

    // matrix
    const T i(0,1) ;
    qclab::dense::SquareMatrix< T >  CY_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 0 ,-i ,
                                               0 , 0 , i , 0 ) ;
    EXPECT_TRUE( cy.matrix() == CY_check ) ;

    // qubit
    EXPECT_EQ( cy.qubit() , 0 ) ;

    // qubits
    auto qubits = cy.qubits() ;
    EXPECT_EQ( qubits.size() , 2 ) ;
    EXPECT_EQ( qubits[0] , 0 ) ;
    EXPECT_EQ( qubits[1] , 1 ) ;
    int qnew[] = { 5 , 3 } ;
    cy.setQubits( &qnew[0] ) ;
    EXPECT_EQ( cy.qubits()[0] , 3 ) ;
    EXPECT_EQ( cy.qubits()[1] , 5 ) ;
    qnew[0] = 0 ;
    qnew[1] = 1 ;
    cy.setQubits( &qnew[0] ) ;

    // print
    cy.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cy.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cy q[0], q[1];\n" ) ;
    std::cout << qasm.str() ;

    // gate
    qclab::qgates::PauliY< T >  Y ;
    EXPECT_TRUE( *cy.gate() == Y ) ;
    EXPECT_TRUE( cy.gate()->matrix() == Y.matrix() ) ;

    // operators == and !=
    EXPECT_TRUE(  cy != Y ) ;
    EXPECT_FALSE( cy == Y ) ;
    qclab::qgates::CY< T >  cy2 ;
    EXPECT_TRUE(  cy == cy2 ) ;
    EXPECT_FALSE( cy != cy2 ) ;

    // setControl, setTarget, setControlState
    cy.setControl( 3 ) ;
    EXPECT_EQ( cy.control() , 3 ) ;
    cy.setTarget( 5 ) ;
    EXPECT_EQ( cy.target() , 5 ) ;
    EXPECT_TRUE(  cy == cy2 ) ;
    EXPECT_FALSE( cy != cy2 ) ;

    cy.setControl( 4 ) ;
    EXPECT_EQ( cy.control() , 4 ) ;
    cy.setTarget( 1 ) ;
    EXPECT_EQ( cy.target() , 1 ) ;
    EXPECT_TRUE(  cy != cy2 ) ;
    EXPECT_FALSE( cy == cy2 ) ;

    qubits[0] = 1 ;
    qubits[1] = 2 ;
    cy.setQubits( &qubits[0] ) ;
    EXPECT_EQ( cy.control() , 1 ) ;
    EXPECT_EQ( cy.target() , 2 ) ;
    EXPECT_TRUE(  cy == cy2 ) ;
    EXPECT_FALSE( cy != cy2 ) ;

    cy.setControl( 0 ) ;
    EXPECT_EQ( cy.control() , 0 ) ;
    cy.setTarget( 1 ) ;
    EXPECT_EQ( cy.target() , 1 ) ;
    cy.setControlState( 0 ) ;
    EXPECT_EQ( cy.controlState() , 0 ) ;
    EXPECT_TRUE(  cy != cy2 ) ;
    EXPECT_FALSE( cy == cy2 ) ;
  }

  {
    qclab::qgates::CY< T >  cy01( 0 , 1 ) ;

    EXPECT_EQ( cy01.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cy01.fixed() ) ;           // fixed
    EXPECT_TRUE( cy01.controlled() ) ;      // controlled
    EXPECT_EQ( cy01.control() , 0 ) ;       // control
    EXPECT_EQ( cy01.target() , 1 ) ;        // target
    EXPECT_EQ( cy01.controlState() , 1 ) ;  // controlState

    // matrix
    const T i(0,1) ;
    qclab::dense::SquareMatrix< T >  CY_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 0 ,-i ,
                                               0 , 0 , i , 0 ) ;
    EXPECT_TRUE( cy01.matrix() == CY_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cy01.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cy q[0], q[1];\n" ) ;
  }

  {
    qclab::qgates::CY< T >  cy35( 3 , 5 ) ;

    EXPECT_EQ( cy35.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cy35.fixed() ) ;           // fixed
    EXPECT_TRUE( cy35.controlled() ) ;      // controlled
    EXPECT_EQ( cy35.control() , 3 ) ;       // control
    EXPECT_EQ( cy35.target() , 5 ) ;        // target
    EXPECT_EQ( cy35.controlState() , 1 ) ;  // controlState

    // matrix
    const T i(0,1) ;
    qclab::dense::SquareMatrix< T >  CY_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 0 ,-i ,
                                               0 , 0 , i , 0 ) ;
    EXPECT_TRUE( cy35.matrix() == CY_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cy35.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cy q[3], q[5];\n" ) ;
  }

  {
    qclab::qgates::CY< T >  cy10( 1 , 0 ) ;

    EXPECT_EQ( cy10.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cy10.fixed() ) ;           // fixed
    EXPECT_TRUE( cy10.controlled() ) ;      // controlled
    EXPECT_EQ( cy10.control() , 1 ) ;       // control
    EXPECT_EQ( cy10.target() , 0 ) ;        // target
    EXPECT_EQ( cy10.controlState() , 1 ) ;  // controlState

    // matrix
    const T i(0,1) ;
    qclab::dense::SquareMatrix< T >  CY_check( 1 , 0 , 0 , 0 ,
                                               0 , 0 , 0 ,-i ,
                                               0 , 0 , 1 , 0 ,
                                               0 , i , 0 , 0 ) ;
    EXPECT_TRUE( cy10.matrix() == CY_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cy10.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cy q[1], q[0];\n" ) ;
  }

  {
    qclab::qgates::CY< T >  cy53( 5 , 3 ) ;

    EXPECT_EQ( cy53.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cy53.fixed() ) ;           // fixed
    EXPECT_TRUE( cy53.controlled() ) ;      // controlled
    EXPECT_EQ( cy53.control() , 5 ) ;       // control
    EXPECT_EQ( cy53.target() , 3 ) ;        // target
    EXPECT_EQ( cy53.controlState() , 1 ) ;  // controlState

    // matrix
    const T i(0,1) ;
    qclab::dense::SquareMatrix< T >  CY_check( 1 , 0 , 0 , 0 ,
                                               0 , 0 , 0 ,-i ,
                                               0 , 0 , 1 , 0 ,
                                               0 , i , 0 , 0 ) ;
    EXPECT_TRUE( cy53.matrix() == CY_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cy53.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cy q[5], q[3];\n" ) ;
  }


  //
  // CYs |0> controlled
  //
  {
    qclab::qgates::CY< T >  cy01( 0 , 1 , 0 ) ;

    EXPECT_EQ( cy01.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cy01.fixed() ) ;           // fixed
    EXPECT_TRUE( cy01.controlled() ) ;      // controlled
    EXPECT_EQ( cy01.control() , 0 ) ;       // control
    EXPECT_EQ( cy01.target() , 1 ) ;        // target
    EXPECT_EQ( cy01.controlState() , 0 ) ;  // controlState

    // matrix
    const T i(0,1) ;
    qclab::dense::SquareMatrix< T >  CY_check( 0 ,-i , 0 , 0 ,
                                               i , 0 , 0 , 0 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 0 , 0 , 1 ) ;
    EXPECT_TRUE( cy01.matrix() == CY_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cy01.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[0];\ncy q[0], q[1];\nx q[0];\n" ) ;
  }

  {
    qclab::qgates::CY< T >  cy10( 1 , 0 , 0 ) ;

    EXPECT_EQ( cy10.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cy10.fixed() ) ;           // fixed
    EXPECT_TRUE( cy10.controlled() ) ;      // controlled
    EXPECT_EQ( cy10.control() , 1 ) ;       // control
    EXPECT_EQ( cy10.target() , 0 ) ;        // target
    EXPECT_EQ( cy10.controlState() , 0 ) ;  // controlState

    // matrix
    const T i(0,1) ;
    qclab::dense::SquareMatrix< T >  CY_check( 0 , 0 ,-i , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               i , 0 , 0 , 0 ,
                                               0 , 0 , 0 , 1 ) ;
    EXPECT_TRUE( cy10.matrix() == CY_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cy10.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[1];\ncy q[1], q[0];\nx q[1];\n" ) ;
  }

  using V = std::vector< T > ;

  const V v2 = { T(1,1) , T(2,2) , T(3,3) , T(4,4) } ;
  const V v3 = { T(1,1) , T(2,2) , T(3,3) , T(4,4) ,
                 T(5,5) , T(6,6) , T(7,7) , T(8,8) } ;
  V vt( 16 ) ; for ( int i = 0; i < 16; ++i ) vt[i] = T( i+1 , i+1 ) ;
  const V v4 = vt ;
  vt = V( 32 ) ; for ( int i = 0; i < 32; ++i ) vt[i] = T( i+1 , i+1 ) ;
  const V v5 = vt ;
  vt = V( 256 ) ; for ( int i = 0; i < 256; ++i ) vt[i] = T( i+1 , i+1 ) ;
  const V v8 = vt ;

  //
  // CNOTs |1> controlled
  //

  // control < target
  {
    qclab::qgates::CY< T >  cy01( 0 , 1 ) ;
    qclab::qgates::CY< T >  cy12( 1 , 2 ) ;
    qclab::qgates::CY< T >  cy02( 0 , 2 ) ;
    qclab::qgates::CY< T >  cy13( 1 , 3 ) ;
    qclab::qgates::CY< T >  cy26( 2 , 6 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cy01.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { T(1,1) , T(2,2) , T(4,-4) , T(-3,3) } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cy01.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cy01.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cy01.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cy01.apply( qclab::Op::Trans , 2 , vec2 ) ;
    check2 = { T(1,1) , T(2,2) , T(-4,4) , T(3,-3) } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cy01.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    cy01.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { T( 1, 1) , T( 2, 2) , T( 3, 3) , T( 4, 4) ,
                 T( 7,-7) , T( 8,-8) , T(-5, 5) , T(-6, 6) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy01.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy01.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy01.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy01.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { T( 1, 1) , T( 2, 2) , T( 3, 3) , T( 4, 4) ,
               T(-7, 7) , T(-8, 8) , T( 5,-5) , T( 6,-6) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy01.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy12.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T( 1, 1) , T( 2, 2) , T( 4,-4) , T(-3, 3) ,
               T( 5, 5) , T( 6, 6) , T( 8,-8) , T(-7, 7) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy12.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy12.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy12.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy12.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { T( 1, 1) , T( 2, 2) , T(-4, 4) , T( 3,-3) ,
               T( 5, 5) , T( 6, 6) , T(-8, 8) , T( 7,-7) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy12.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy02.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T( 1, 1) , T( 2, 2) , T( 3, 3) , T( 4, 4) ,
               T( 6,-6) , T(-5, 5) , T( 8,-8) , T(-7, 7) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy02.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy02.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy02.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy02.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { T( 1, 1) , T( 2, 2) , T( 3, 3) , T( 4, 4) ,
               T(-6, 6) , T( 5,-5) , T(-8, 8) , T( 7,-7) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy02.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    cy12.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { T(  1,  1) , T(  2,  2) , T(  3,  3) , T(  4,  4) ,
                 T(  7, -7) , T(  8, -8) , T( -5,  5) , T( -6,  6) ,
                 T(  9,  9) , T( 10, 10) , T( 11, 11) , T( 12, 12) ,
                 T( 15,-15) , T( 16,-16) , T(-13, 13) , T(-14, 14) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cy12.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cy12.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cy12.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cy12.apply( qclab::Op::Trans , 4 , vec4 ) ;
    check4 = { T(  1,  1) , T(  2,  2) , T(  3,  3) , T(  4,  4) ,
               T( -7,  7) , T( -8,  8) , T(  5, -5) , T(  6, -6) ,
               T(  9,  9) , T( 10, 10) , T( 11, 11) , T( 12, 12) ,
               T(-15, 15) , T(-16, 16) , T( 13,-13) , T( 14,-14) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cy12.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    cy13.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = { T(  1,  1) , T(  2,  2) , T(  3,  3) , T(  4,  4) ,
                 T(  5,  5) , T(  6,  6) , T(  7,  7) , T(  8,  8) ,
                 T( 11,-11) , T( 12,-12) , T( -9,  9) , T(-10, 10) ,
                 T( 15,-15) , T( 16,-16) , T(-13, 13) , T(-14, 14) ,
                 T( 17, 17) , T( 18, 18) , T( 19, 19) , T( 20, 20) ,
                 T( 21, 21) , T( 22, 22) , T( 23, 23) , T( 24, 24) ,
                 T( 27,-27) , T( 28,-28) , T(-25, 25) , T(-26, 26) ,
                 T( 31,-31) , T( 32,-32) , T(-29, 29) , T(-30, 30) } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cy13.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cy13.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cy13.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cy13.apply( qclab::Op::Trans , 5 , vec5 ) ;
    check5 = { T(  1,  1) , T(  2,  2) , T(  3,  3) , T(  4,  4) ,
               T(  5,  5) , T(  6,  6) , T(  7,  7) , T(  8,  8) ,
               T(-11, 11) , T(-12, 12) , T(  9, -9) , T( 10,-10) ,
               T(-15, 15) , T(-16, 16) , T( 13,-13) , T( 14,-14) ,
               T( 17, 17) , T( 18, 18) , T( 19, 19) , T( 20, 20) ,
               T( 21, 21) , T( 22, 22) , T( 23, 23) , T( 24, 24) ,
               T(-27, 27) , T(-28, 28) , T( 25,-25) , T( 26,-26) ,
               T(-31, 31) , T(-32, 32) , T( 29,-29) , T( 30,-30) } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cy13.apply_device( qclab::Op::Trans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    cy26.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    V check8 = { T(   1,   1) , T(   2,   2) , T(   3,   3) , T(   4,   4) ,
                 T(   5,   5) , T(   6,   6) , T(   7,   7) , T(   8,   8) ,
                 T(   9,   9) , T(  10,  10) , T(  11,  11) , T(  12,  12) ,
                 T(  13,  13) , T(  14,  14) , T(  15,  15) , T(  16,  16) ,
                 T(  17,  17) , T(  18,  18) , T(  19,  19) , T(  20,  20) ,
                 T(  21,  21) , T(  22,  22) , T(  23,  23) , T(  24,  24) ,
                 T(  25,  25) , T(  26,  26) , T(  27,  27) , T(  28,  28) ,
                 T(  29,  29) , T(  30,  30) , T(  31,  31) , T(  32,  32) ,
                 T(  35, -35) , T(  36, -36) , T( -33,  33) , T( -34,  34) ,
                 T(  39, -39) , T(  40, -40) , T( -37,  37) , T( -38,  38) ,
                 T(  43, -43) , T(  44, -44) , T( -41,  41) , T( -42,  42) ,
                 T(  47, -47) , T(  48, -48) , T( -45,  45) , T( -46,  46) ,
                 T(  51, -51) , T(  52, -52) , T( -49,  49) , T( -50,  50) ,
                 T(  55, -55) , T(  56, -56) , T( -53,  53) , T( -54,  54) ,
                 T(  59, -59) , T(  60, -60) , T( -57,  57) , T( -58,  58) ,
                 T(  63, -63) , T(  64, -64) , T( -61,  61) , T( -62,  62) ,
                 T(  65,  65) , T(  66,  66) , T(  67,  67) , T(  68,  68) ,
                 T(  69,  69) , T(  70,  70) , T(  71,  71) , T(  72,  72) ,
                 T(  73,  73) , T(  74,  74) , T(  75,  75) , T(  76,  76) ,
                 T(  77,  77) , T(  78,  78) , T(  79,  79) , T(  80,  80) ,
                 T(  81,  81) , T(  82,  82) , T(  83,  83) , T(  84,  84) ,
                 T(  85,  85) , T(  86,  86) , T(  87,  87) , T(  88,  88) ,
                 T(  89,  89) , T(  90,  90) , T(  91,  91) , T(  92,  92) ,
                 T(  93,  93) , T(  94,  94) , T(  95,  95) , T(  96,  96) ,
                 T(  99, -99) , T( 100,-100) , T( -97,  97) , T( -98,  98) ,
                 T( 103,-103) , T( 104,-104) , T(-101, 101) , T(-102, 102) ,
                 T( 107,-107) , T( 108,-108) , T(-105, 105) , T(-106, 106) ,
                 T( 111,-111) , T( 112,-112) , T(-109, 109) , T(-110, 110) ,
                 T( 115,-115) , T( 116,-116) , T(-113, 113) , T(-114, 114) ,
                 T( 119,-119) , T( 120,-120) , T(-117, 117) , T(-118, 118) ,
                 T( 123,-123) , T( 124,-124) , T(-121, 121) , T(-122, 122) ,
                 T( 127,-127) , T( 128,-128) , T(-125, 125) , T(-126, 126) ,
                 T( 129, 129) , T( 130, 130) , T( 131, 131) , T( 132, 132) ,
                 T( 133, 133) , T( 134, 134) , T( 135, 135) , T( 136, 136) ,
                 T( 137, 137) , T( 138, 138) , T( 139, 139) , T( 140, 140) ,
                 T( 141, 141) , T( 142, 142) , T( 143, 143) , T( 144, 144) ,
                 T( 145, 145) , T( 146, 146) , T( 147, 147) , T( 148, 148) ,
                 T( 149, 149) , T( 150, 150) , T( 151, 151) , T( 152, 152) ,
                 T( 153, 153) , T( 154, 154) , T( 155, 155) , T( 156, 156) ,
                 T( 157, 157) , T( 158, 158) , T( 159, 159) , T( 160, 160) ,
                 T( 163,-163) , T( 164,-164) , T(-161, 161) , T(-162, 162) ,
                 T( 167,-167) , T( 168,-168) , T(-165, 165) , T(-166, 166) ,
                 T( 171,-171) , T( 172,-172) , T(-169, 169) , T(-170, 170) ,
                 T( 175,-175) , T( 176,-176) , T(-173, 173) , T(-174, 174) ,
                 T( 179,-179) , T( 180,-180) , T(-177, 177) , T(-178, 178) ,
                 T( 183,-183) , T( 184,-184) , T(-181, 181) , T(-182, 182) ,
                 T( 187,-187) , T( 188,-188) , T(-185, 185) , T(-186, 186) ,
                 T( 191,-191) , T( 192,-192) , T(-189, 189) , T(-190, 190) ,
                 T( 193, 193) , T( 194, 194) , T( 195, 195) , T( 196, 196) ,
                 T( 197, 197) , T( 198, 198) , T( 199, 199) , T( 200, 200) ,
                 T( 201, 201) , T( 202, 202) , T( 203, 203) , T( 204, 204) ,
                 T( 205, 205) , T( 206, 206) , T( 207, 207) , T( 208, 208) ,
                 T( 209, 209) , T( 210, 210) , T( 211, 211) , T( 212, 212) ,
                 T( 213, 213) , T( 214, 214) , T( 215, 215) , T( 216, 216) ,
                 T( 217, 217) , T( 218, 218) , T( 219, 219) , T( 220, 220) ,
                 T( 221, 221) , T( 222, 222) , T( 223, 223) , T( 224, 224) ,
                 T( 227,-227) , T( 228,-228) , T(-225, 225) , T(-226, 226) ,
                 T( 231,-231) , T( 232,-232) , T(-229, 229) , T(-230, 230) ,
                 T( 235,-235) , T( 236,-236) , T(-233, 233) , T(-234, 234) ,
                 T( 239,-239) , T( 240,-240) , T(-237, 237) , T(-238, 238) ,
                 T( 243,-243) , T( 244,-244) , T(-241, 241) , T(-242, 242) ,
                 T( 247,-247) , T( 248,-248) , T(-245, 245) , T(-246, 246) ,
                 T( 251,-251) , T( 252,-252) , T(-249, 249) , T(-250, 250) ,
                 T( 255,-255) , T( 256,-256) , T(-253, 253) , T(-254, 254) } ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cy26.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cy26.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cy26.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cy26.apply( qclab::Op::Trans , 8 , vec8 ) ;
    check8 = { T(   1,   1) , T(   2,   2) , T(   3,   3) , T(   4,   4) ,
               T(   5,   5) , T(   6,   6) , T(   7,   7) , T(   8,   8) ,
               T(   9,   9) , T(  10,  10) , T(  11,  11) , T(  12,  12) ,
               T(  13,  13) , T(  14,  14) , T(  15,  15) , T(  16,  16) ,
               T(  17,  17) , T(  18,  18) , T(  19,  19) , T(  20,  20) ,
               T(  21,  21) , T(  22,  22) , T(  23,  23) , T(  24,  24) ,
               T(  25,  25) , T(  26,  26) , T(  27,  27) , T(  28,  28) ,
               T(  29,  29) , T(  30,  30) , T(  31,  31) , T(  32,  32) ,
               T( -35,  35) , T( -36,  36) , T(  33, -33) , T(  34, -34) ,
               T( -39,  39) , T( -40,  40) , T(  37, -37) , T(  38, -38) ,
               T( -43,  43) , T( -44,  44) , T(  41, -41) , T(  42, -42) ,
               T( -47,  47) , T( -48,  48) , T(  45, -45) , T(  46, -46) ,
               T( -51,  51) , T( -52,  52) , T(  49, -49) , T(  50, -50) ,
               T( -55,  55) , T( -56,  56) , T(  53, -53) , T(  54, -54) ,
               T( -59,  59) , T( -60,  60) , T(  57, -57) , T(  58, -58) ,
               T( -63,  63) , T( -64,  64) , T(  61, -61) , T(  62, -62) ,
               T(  65,  65) , T(  66,  66) , T(  67,  67) , T(  68,  68) ,
               T(  69,  69) , T(  70,  70) , T(  71,  71) , T(  72,  72) ,
               T(  73,  73) , T(  74,  74) , T(  75,  75) , T(  76,  76) ,
               T(  77,  77) , T(  78,  78) , T(  79,  79) , T(  80,  80) ,
               T(  81,  81) , T(  82,  82) , T(  83,  83) , T(  84,  84) ,
               T(  85,  85) , T(  86,  86) , T(  87,  87) , T(  88,  88) ,
               T(  89,  89) , T(  90,  90) , T(  91,  91) , T(  92,  92) ,
               T(  93,  93) , T(  94,  94) , T(  95,  95) , T(  96,  96) ,
               T( -99,  99) , T(-100, 100) , T(  97, -97) , T(  98, -98) ,
               T(-103, 103) , T(-104, 104) , T( 101,-101) , T( 102,-102) ,
               T(-107, 107) , T(-108, 108) , T( 105,-105) , T( 106,-106) ,
               T(-111, 111) , T(-112, 112) , T( 109,-109) , T( 110,-110) ,
               T(-115, 115) , T(-116, 116) , T( 113,-113) , T( 114,-114) ,
               T(-119, 119) , T(-120, 120) , T( 117,-117) , T( 118,-118) ,
               T(-123, 123) , T(-124, 124) , T( 121,-121) , T( 122,-122) ,
               T(-127, 127) , T(-128, 128) , T( 125,-125) , T( 126,-126) ,
               T( 129, 129) , T( 130, 130) , T( 131, 131) , T( 132, 132) ,
               T( 133, 133) , T( 134, 134) , T( 135, 135) , T( 136, 136) ,
               T( 137, 137) , T( 138, 138) , T( 139, 139) , T( 140, 140) ,
               T( 141, 141) , T( 142, 142) , T( 143, 143) , T( 144, 144) ,
               T( 145, 145) , T( 146, 146) , T( 147, 147) , T( 148, 148) ,
               T( 149, 149) , T( 150, 150) , T( 151, 151) , T( 152, 152) ,
               T( 153, 153) , T( 154, 154) , T( 155, 155) , T( 156, 156) ,
               T( 157, 157) , T( 158, 158) , T( 159, 159) , T( 160, 160) ,
               T(-163, 163) , T(-164, 164) , T( 161,-161) , T( 162,-162) ,
               T(-167, 167) , T(-168, 168) , T( 165,-165) , T( 166,-166) ,
               T(-171, 171) , T(-172, 172) , T( 169,-169) , T( 170,-170) ,
               T(-175, 175) , T(-176, 176) , T( 173,-173) , T( 174,-174) ,
               T(-179, 179) , T(-180, 180) , T( 177,-177) , T( 178,-178) ,
               T(-183, 183) , T(-184, 184) , T( 181,-181) , T( 182,-182) ,
               T(-187, 187) , T(-188, 188) , T( 185,-185) , T( 186,-186) ,
               T(-191, 191) , T(-192, 192) , T( 189,-189) , T( 190,-190) ,
               T( 193, 193) , T( 194, 194) , T( 195, 195) , T( 196, 196) ,
               T( 197, 197) , T( 198, 198) , T( 199, 199) , T( 200, 200) ,
               T( 201, 201) , T( 202, 202) , T( 203, 203) , T( 204, 204) ,
               T( 205, 205) , T( 206, 206) , T( 207, 207) , T( 208, 208) ,
               T( 209, 209) , T( 210, 210) , T( 211, 211) , T( 212, 212) ,
               T( 213, 213) , T( 214, 214) , T( 215, 215) , T( 216, 216) ,
               T( 217, 217) , T( 218, 218) , T( 219, 219) , T( 220, 220) ,
               T( 221, 221) , T( 222, 222) , T( 223, 223) , T( 224, 224) ,
               T(-227, 227) , T(-228, 228) , T( 225,-225) , T( 226,-226) ,
               T(-231, 231) , T(-232, 232) , T( 229,-229) , T( 230,-230) ,
               T(-235, 235) , T(-236, 236) , T( 233,-233) , T( 234,-234) ,
               T(-239, 239) , T(-240, 240) , T( 237,-237) , T( 238,-238) ,
               T(-243, 243) , T(-244, 244) , T( 241,-241) , T( 242,-242) ,
               T(-247, 247) , T(-248, 248) , T( 245,-245) , T( 246,-246) ,
               T(-251, 251) , T(-252, 252) , T( 249,-249) , T( 250,-250) ,
               T(-255, 255) , T(-256, 256) , T( 253,-253) , T( 254,-254) } ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cy26.apply_device( qclab::Op::Trans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

  // control > target
  {
    qclab::qgates::CY< T >  cy10( 1 , 0 ) ;
    qclab::qgates::CY< T >  cy21( 2 , 1 ) ;
    qclab::qgates::CY< T >  cy20( 2 , 0 ) ;
    qclab::qgates::CY< T >  cy31( 3 , 1 ) ;
    qclab::qgates::CY< T >  cy62( 6 , 2 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cy10.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { T(1,1) , T(4,-4) , T(3,3) , T(-2,2) } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cy10.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cy10.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cy10.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cy10.apply( qclab::Op::Trans , 2 , vec2 ) ;
    check2 = { T(1,1) , T(-4,4) , T(3,3) , T(2,-2) } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cy10.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    cy10.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { T( 1, 1) , T( 2, 2) , T( 7,-7) , T( 8,-8) ,
                 T( 5, 5) , T( 6, 6) , T(-3, 3) , T(-4, 4) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy10.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy10.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy10.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy10.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { T( 1, 1) , T( 2, 2) , T(-7, 7) , T(-8, 8) ,
               T( 5, 5) , T( 6, 6) , T( 3,-3) , T( 4,-4) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy10.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy21.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T( 1, 1) , T( 4,-4) , T( 3, 3) , T(-2, 2) ,
               T( 5, 5) , T( 8,-8) , T( 7, 7) , T(-6, 6) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy21.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy21.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy21.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy21.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { T( 1, 1) , T(-4, 4) , T( 3, 3) , T( 2,-2) ,
               T( 5, 5) , T(-8, 8) , T( 7, 7) , T( 6,-6) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy21.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy20.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T( 1, 1) , T( 6,-6) , T( 3, 3) , T( 8,-8) ,
               T( 5, 5) , T(-2, 2) , T( 7, 7) , T(-4, 4) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy20.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy20.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy20.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy20.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { T( 1, 1) , T(-6, 6) , T( 3, 3) , T(-8, 8) ,
               T( 5, 5) , T( 2,-2) , T( 7, 7) , T( 4,-4) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy20.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    cy21.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { T(  1,  1) , T(  2,  2) , T(  7, -7) , T(  8, -8) ,
                 T(  5,  5) , T(  6,  6) , T( -3,  3) , T( -4,  4) ,
                 T(  9,  9) , T( 10, 10) , T( 15,-15) , T( 16,-16) ,
                 T( 13, 13) , T( 14, 14) , T(-11, 11) , T(-12, 12) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cy21.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cy21.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cy21.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cy21.apply( qclab::Op::Trans , 4 , vec4 ) ;
    check4 = { T(  1,  1) , T(  2,  2) , T( -7,  7) , T( -8,  8) ,
               T(  5,  5) , T(  6,  6) , T(  3, -3) , T(  4, -4) ,
               T(  9,  9) , T( 10, 10) , T(-15, 15) , T(-16, 16) ,
               T( 13, 13) , T( 14, 14) , T( 11,-11) , T( 12,-12) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cy21.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    cy31.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = { T(  1,  1) , T(  2,  2) , T( 11,-11) , T( 12,-12) ,
                 T(  5,  5) , T(  6,  6) , T( 15,-15) , T( 16,-16) ,
                 T(  9,  9) , T( 10, 10) , T( -3,  3) , T( -4,  4) ,
                 T( 13, 13) , T( 14, 14) , T( -7,  7) , T( -8,  8) ,
                 T( 17, 17) , T( 18, 18) , T( 27,-27) , T( 28,-28) ,
                 T( 21, 21) , T( 22, 22) , T( 31,-31) , T( 32,-32) ,
                 T( 25, 25) , T( 26, 26) , T(-19, 19) , T(-20, 20) ,
                 T( 29, 29) , T( 30, 30) , T(-23, 23) , T(-24, 24) } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cy31.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cy31.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cy31.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cy31.apply( qclab::Op::Trans , 5 , vec5 ) ;
    check5 = { T(  1,  1) , T(  2,  2) , T(-11, 11) , T(-12, 12) ,
               T(  5,  5) , T(  6,  6) , T(-15, 15) , T(-16, 16) ,
               T(  9,  9) , T( 10, 10) , T(  3, -3) , T(  4, -4) ,
               T( 13, 13) , T( 14, 14) , T(  7, -7) , T(  8, -8) ,
               T( 17, 17) , T( 18, 18) , T(-27, 27) , T(-28, 28) ,
               T( 21, 21) , T( 22, 22) , T(-31, 31) , T(-32, 32) ,
               T( 25, 25) , T( 26, 26) , T( 19,-19) , T( 20,-20) ,
               T( 29, 29) , T( 30, 30) , T( 23,-23) , T( 24,-24) } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cy31.apply_device( qclab::Op::Trans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    cy62.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    V check8 = { T(   1,   1) , T(   2,   2) , T(  35, -35) , T(  36, -36) ,
                 T(   5,   5) , T(   6,   6) , T(  39, -39) , T(  40, -40) ,
                 T(   9,   9) , T(  10,  10) , T(  43, -43) , T(  44, -44) ,
                 T(  13,  13) , T(  14,  14) , T(  47, -47) , T(  48, -48) ,
                 T(  17,  17) , T(  18,  18) , T(  51, -51) , T(  52, -52) ,
                 T(  21,  21) , T(  22,  22) , T(  55, -55) , T(  56, -56) ,
                 T(  25,  25) , T(  26,  26) , T(  59, -59) , T(  60, -60) ,
                 T(  29,  29) , T(  30,  30) , T(  63, -63) , T(  64, -64) ,
                 T(  33,  33) , T(  34,  34) , T(  -3,   3) , T(  -4,   4) ,
                 T(  37,  37) , T(  38,  38) , T(  -7,   7) , T(  -8,   8) ,
                 T(  41,  41) , T(  42,  42) , T( -11,  11) , T( -12,  12) ,
                 T(  45,  45) , T(  46,  46) , T( -15,  15) , T( -16,  16) ,
                 T(  49,  49) , T(  50,  50) , T( -19,  19) , T( -20,  20) ,
                 T(  53,  53) , T(  54,  54) , T( -23,  23) , T( -24,  24) ,
                 T(  57,  57) , T(  58,  58) , T( -27,  27) , T( -28,  28) ,
                 T(  61,  61) , T(  62,  62) , T( -31,  31) , T( -32,  32) ,
                 T(  65,  65) , T(  66,  66) , T(  99, -99) , T( 100,-100) ,
                 T(  69,  69) , T(  70,  70) , T( 103,-103) , T( 104,-104) ,
                 T(  73,  73) , T(  74,  74) , T( 107,-107) , T( 108,-108) ,
                 T(  77,  77) , T(  78,  78) , T( 111,-111) , T( 112,-112) ,
                 T(  81,  81) , T(  82,  82) , T( 115,-115) , T( 116,-116) ,
                 T(  85,  85) , T(  86,  86) , T( 119,-119) , T( 120,-120) ,
                 T(  89,  89) , T(  90,  90) , T( 123,-123) , T( 124,-124) ,
                 T(  93,  93) , T(  94,  94) , T( 127,-127) , T( 128,-128) ,
                 T(  97,  97) , T(  98,  98) , T( -67,  67) , T( -68,  68) ,
                 T( 101, 101) , T( 102, 102) , T( -71,  71) , T( -72,  72) ,
                 T( 105, 105) , T( 106, 106) , T( -75,  75) , T( -76,  76) ,
                 T( 109, 109) , T( 110, 110) , T( -79,  79) , T( -80,  80) ,
                 T( 113, 113) , T( 114, 114) , T( -83,  83) , T( -84,  84) ,
                 T( 117, 117) , T( 118, 118) , T( -87,  87) , T( -88,  88) ,
                 T( 121, 121) , T( 122, 122) , T( -91,  91) , T( -92,  92) ,
                 T( 125, 125) , T( 126, 126) , T( -95,  95) , T( -96,  96) ,
                 T( 129, 129) , T( 130, 130) , T( 163,-163) , T( 164,-164) ,
                 T( 133, 133) , T( 134, 134) , T( 167,-167) , T( 168,-168) ,
                 T( 137, 137) , T( 138, 138) , T( 171,-171) , T( 172,-172) ,
                 T( 141, 141) , T( 142, 142) , T( 175,-175) , T( 176,-176) ,
                 T( 145, 145) , T( 146, 146) , T( 179,-179) , T( 180,-180) ,
                 T( 149, 149) , T( 150, 150) , T( 183,-183) , T( 184,-184) ,
                 T( 153, 153) , T( 154, 154) , T( 187,-187) , T( 188,-188) ,
                 T( 157, 157) , T( 158, 158) , T( 191,-191) , T( 192,-192) ,
                 T( 161, 161) , T( 162, 162) , T(-131, 131) , T(-132, 132) ,
                 T( 165, 165) , T( 166, 166) , T(-135, 135) , T(-136, 136) ,
                 T( 169, 169) , T( 170, 170) , T(-139, 139) , T(-140, 140) ,
                 T( 173, 173) , T( 174, 174) , T(-143, 143) , T(-144, 144) ,
                 T( 177, 177) , T( 178, 178) , T(-147, 147) , T(-148, 148) ,
                 T( 181, 181) , T( 182, 182) , T(-151, 151) , T(-152, 152) ,
                 T( 185, 185) , T( 186, 186) , T(-155, 155) , T(-156, 156) ,
                 T( 189, 189) , T( 190, 190) , T(-159, 159) , T(-160, 160) ,
                 T( 193, 193) , T( 194, 194) , T( 227,-227) , T( 228,-228) ,
                 T( 197, 197) , T( 198, 198) , T( 231,-231) , T( 232,-232) ,
                 T( 201, 201) , T( 202, 202) , T( 235,-235) , T( 236,-236) ,
                 T( 205, 205) , T( 206, 206) , T( 239,-239) , T( 240,-240) ,
                 T( 209, 209) , T( 210, 210) , T( 243,-243) , T( 244,-244) ,
                 T( 213, 213) , T( 214, 214) , T( 247,-247) , T( 248,-248) ,
                 T( 217, 217) , T( 218, 218) , T( 251,-251) , T( 252,-252) ,
                 T( 221, 221) , T( 222, 222) , T( 255,-255) , T( 256,-256) ,
                 T( 225, 225) , T( 226, 226) , T(-195, 195) , T(-196, 196) ,
                 T( 229, 229) , T( 230, 230) , T(-199, 199) , T(-200, 200) ,
                 T( 233, 233) , T( 234, 234) , T(-203, 203) , T(-204, 204) ,
                 T( 237, 237) , T( 238, 238) , T(-207, 207) , T(-208, 208) ,
                 T( 241, 241) , T( 242, 242) , T(-211, 211) , T(-212, 212) ,
                 T( 245, 245) , T( 246, 246) , T(-215, 215) , T(-216, 216) ,
                 T( 249, 249) , T( 250, 250) , T(-219, 219) , T(-220, 220) ,
                 T( 253, 253) , T( 254, 254) , T(-223, 223) , T(-224, 224) } ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cy62.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cy62.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cy62.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cy62.apply( qclab::Op::Trans , 8 , vec8 ) ;
    check8 = { T(   1,   1) , T(   2,   2) , T( -35,  35) , T( -36,  36) ,
               T(   5,   5) , T(   6,   6) , T( -39,  39) , T( -40,  40) ,
               T(   9,   9) , T(  10,  10) , T( -43,  43) , T( -44,  44) ,
               T(  13,  13) , T(  14,  14) , T( -47,  47) , T( -48,  48) ,
               T(  17,  17) , T(  18,  18) , T( -51,  51) , T( -52,  52) ,
               T(  21,  21) , T(  22,  22) , T( -55,  55) , T( -56,  56) ,
               T(  25,  25) , T(  26,  26) , T( -59,  59) , T( -60,  60) ,
               T(  29,  29) , T(  30,  30) , T( -63,  63) , T( -64,  64) ,
               T(  33,  33) , T(  34,  34) , T(   3,  -3) , T(   4,  -4) ,
               T(  37,  37) , T(  38,  38) , T(   7,  -7) , T(   8,  -8) ,
               T(  41,  41) , T(  42,  42) , T(  11, -11) , T(  12, -12) ,
               T(  45,  45) , T(  46,  46) , T(  15, -15) , T(  16, -16) ,
               T(  49,  49) , T(  50,  50) , T(  19, -19) , T(  20, -20) ,
               T(  53,  53) , T(  54,  54) , T(  23, -23) , T(  24, -24) ,
               T(  57,  57) , T(  58,  58) , T(  27, -27) , T(  28, -28) ,
               T(  61,  61) , T(  62,  62) , T(  31, -31) , T(  32, -32) ,
               T(  65,  65) , T(  66,  66) , T( -99,  99) , T(-100, 100) ,
               T(  69,  69) , T(  70,  70) , T(-103, 103) , T(-104, 104) ,
               T(  73,  73) , T(  74,  74) , T(-107, 107) , T(-108, 108) ,
               T(  77,  77) , T(  78,  78) , T(-111, 111) , T(-112, 112) ,
               T(  81,  81) , T(  82,  82) , T(-115, 115) , T(-116, 116) ,
               T(  85,  85) , T(  86,  86) , T(-119, 119) , T(-120, 120) ,
               T(  89,  89) , T(  90,  90) , T(-123, 123) , T(-124, 124) ,
               T(  93,  93) , T(  94,  94) , T(-127, 127) , T(-128, 128) ,
               T(  97,  97) , T(  98,  98) , T(  67, -67) , T(  68, -68) ,
               T( 101, 101) , T( 102, 102) , T(  71, -71) , T(  72, -72) ,
               T( 105, 105) , T( 106, 106) , T(  75, -75) , T(  76, -76) ,
               T( 109, 109) , T( 110, 110) , T(  79, -79) , T(  80, -80) ,
               T( 113, 113) , T( 114, 114) , T(  83, -83) , T(  84, -84) ,
               T( 117, 117) , T( 118, 118) , T(  87, -87) , T(  88, -88) ,
               T( 121, 121) , T( 122, 122) , T(  91, -91) , T(  92, -92) ,
               T( 125, 125) , T( 126, 126) , T(  95, -95) , T(  96, -96) ,
               T( 129, 129) , T( 130, 130) , T(-163, 163) , T(-164, 164) ,
               T( 133, 133) , T( 134, 134) , T(-167, 167) , T(-168, 168) ,
               T( 137, 137) , T( 138, 138) , T(-171, 171) , T(-172, 172) ,
               T( 141, 141) , T( 142, 142) , T(-175, 175) , T(-176, 176) ,
               T( 145, 145) , T( 146, 146) , T(-179, 179) , T(-180, 180) ,
               T( 149, 149) , T( 150, 150) , T(-183, 183) , T(-184, 184) ,
               T( 153, 153) , T( 154, 154) , T(-187, 187) , T(-188, 188) ,
               T( 157, 157) , T( 158, 158) , T(-191, 191) , T(-192, 192) ,
               T( 161, 161) , T( 162, 162) , T( 131,-131) , T( 132,-132) ,
               T( 165, 165) , T( 166, 166) , T( 135,-135) , T( 136,-136) ,
               T( 169, 169) , T( 170, 170) , T( 139,-139) , T( 140,-140) ,
               T( 173, 173) , T( 174, 174) , T( 143,-143) , T( 144,-144) ,
               T( 177, 177) , T( 178, 178) , T( 147,-147) , T( 148,-148) ,
               T( 181, 181) , T( 182, 182) , T( 151,-151) , T( 152,-152) ,
               T( 185, 185) , T( 186, 186) , T( 155,-155) , T( 156,-156) ,
               T( 189, 189) , T( 190, 190) , T( 159,-159) , T( 160,-160) ,
               T( 193, 193) , T( 194, 194) , T(-227, 227) , T(-228, 228) ,
               T( 197, 197) , T( 198, 198) , T(-231, 231) , T(-232, 232) ,
               T( 201, 201) , T( 202, 202) , T(-235, 235) , T(-236, 236) ,
               T( 205, 205) , T( 206, 206) , T(-239, 239) , T(-240, 240) ,
               T( 209, 209) , T( 210, 210) , T(-243, 243) , T(-244, 244) ,
               T( 213, 213) , T( 214, 214) , T(-247, 247) , T(-248, 248) ,
               T( 217, 217) , T( 218, 218) , T(-251, 251) , T(-252, 252) ,
               T( 221, 221) , T( 222, 222) , T(-255, 255) , T(-256, 256) ,
               T( 225, 225) , T( 226, 226) , T( 195,-195) , T( 196,-196) ,
               T( 229, 229) , T( 230, 230) , T( 199,-199) , T( 200,-200) ,
               T( 233, 233) , T( 234, 234) , T( 203,-203) , T( 204,-204) ,
               T( 237, 237) , T( 238, 238) , T( 207,-207) , T( 208,-208) ,
               T( 241, 241) , T( 242, 242) , T( 211,-211) , T( 212,-212) ,
               T( 245, 245) , T( 246, 246) , T( 215,-215) , T( 216,-216) ,
               T( 249, 249) , T( 250, 250) , T( 219,-219) , T( 220,-220) ,
               T( 253, 253) , T( 254, 254) , T( 223,-223) , T( 224,-224) } ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cy62.apply_device( qclab::Op::Trans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

  //
  // CNOTs |0> controlled
  //

  // control < target
  {
    qclab::qgates::CY< T >  cy01( 0 , 1 , 0 ) ;
    qclab::qgates::CY< T >  cy12( 1 , 2 , 0 ) ;
    qclab::qgates::CY< T >  cy02( 0 , 2 , 0 ) ;
    qclab::qgates::CY< T >  cy13( 1 , 3 , 0 ) ;
    qclab::qgates::CY< T >  cy26( 2 , 6 , 0 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cy01.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { T(2,-2) , T(-1,1) , T(3,3) , T(4,4) } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cy01.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cy01.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cy01.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cy01.apply( qclab::Op::Trans , 2 , vec2 ) ;
    check2 = { T(-2,2) , T(1,-1) , T(3,3) , T(4,4) } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cy01.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    cy01.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { T( 3,-3) , T( 4,-4) , T(-1, 1) , T(-2, 2) ,
                 T( 5, 5) , T( 6, 6) , T( 7, 7) , T( 8, 8) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy01.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy01.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy01.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy01.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { T(-3, 3) , T(-4, 4) , T( 1,-1) , T( 2,-2) ,
               T( 5, 5) , T( 6, 6) , T( 7, 7) , T( 8, 8) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy01.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy12.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T( 2,-2) , T(-1, 1) , T( 3, 3) , T( 4, 4) ,
               T( 6,-6) , T(-5, 5) , T( 7, 7) , T( 8, 8) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy12.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy12.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy12.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy12.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { T(-2, 2) , T( 1,-1) , T( 3, 3) , T( 4, 4) ,
               T(-6, 6) , T( 5,-5) , T( 7, 7) , T( 8, 8) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy12.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy02.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T( 2,-2) , T(-1, 1) , T( 4,-4) , T(-3, 3) ,
               T( 5, 5) , T( 6, 6) , T( 7, 7) , T( 8, 8) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy02.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy02.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy02.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy02.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { T(-2, 2) , T( 1,-1) , T(-4, 4) , T( 3,-3) ,
               T( 5, 5) , T( 6, 6) , T( 7, 7) , T( 8, 8) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy02.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    cy12.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { T(  3, -3) , T(  4, -4) , T( -1,  1) , T( -2,  2) ,
                 T(  5,  5) , T(  6,  6) , T(  7,  7) , T(  8,  8) ,
                 T( 11,-11) , T( 12,-12) , T( -9,  9) , T(-10, 10) ,
                 T( 13, 13) , T( 14, 14) , T( 15, 15) , T( 16, 16) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cy12.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cy12.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cy12.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cy12.apply( qclab::Op::Trans , 4 , vec4 ) ;
    check4 = { T( -3,  3) , T( -4,  4) , T(  1, -1) , T(  2, -2) ,
               T(  5,  5) , T(  6,  6) , T(  7,  7) , T(  8,  8) ,
               T(-11, 11) , T(-12, 12) , T(  9, -9) , T( 10,-10) ,
               T( 13, 13) , T( 14, 14) , T( 15, 15) , T( 16, 16) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cy12.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    cy13.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = { T(  3, -3) , T(  4, -4) , T( -1,  1) , T( -2,  2) ,
                 T(  7, -7) , T(  8, -8) , T( -5,  5) , T( -6,  6) ,
                 T(  9,  9) , T( 10, 10) , T( 11, 11) , T( 12, 12) ,
                 T( 13, 13) , T( 14, 14) , T( 15, 15) , T( 16, 16) ,
                 T( 19,-19) , T( 20,-20) , T(-17, 17) , T(-18, 18) ,
                 T( 23,-23) , T( 24,-24) , T(-21, 21) , T(-22, 22) ,
                 T( 25, 25) , T( 26, 26) , T( 27, 27) , T( 28, 28) ,
                 T( 29, 29) , T( 30, 30) , T( 31, 31) , T( 32, 32) } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cy13.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cy13.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cy13.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cy13.apply( qclab::Op::Trans , 5 , vec5 ) ;
    check5 = { T( -3,  3) , T( -4,  4) , T(  1, -1) , T(  2, -2) ,
               T( -7,  7) , T( -8,  8) , T(  5, -5) , T(  6, -6) ,
               T(  9,  9) , T( 10, 10) , T( 11, 11) , T( 12, 12) ,
               T( 13, 13) , T( 14, 14) , T( 15, 15) , T( 16, 16) ,
               T(-19, 19) , T(-20, 20) , T( 17,-17) , T( 18,-18) ,
               T(-23, 23) , T(-24, 24) , T( 21,-21) , T( 22,-22) ,
               T( 25, 25) , T( 26, 26) , T( 27, 27) , T( 28, 28) ,
               T( 29, 29) , T( 30, 30) , T( 31, 31) , T( 32, 32) } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cy13.apply_device( qclab::Op::Trans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    cy26.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    V check8 = { T(   3,  -3) , T(   4,  -4) , T(  -1,   1) , T(  -2,   2) ,
                 T(   7,  -7) , T(   8,  -8) , T(  -5,   5) , T(  -6,   6) ,
                 T(  11, -11) , T(  12, -12) , T(  -9,   9) , T( -10,  10) ,
                 T(  15, -15) , T(  16, -16) , T( -13,  13) , T( -14,  14) ,
                 T(  19, -19) , T(  20, -20) , T( -17,  17) , T( -18,  18) ,
                 T(  23, -23) , T(  24, -24) , T( -21,  21) , T( -22,  22) ,
                 T(  27, -27) , T(  28, -28) , T( -25,  25) , T( -26,  26) ,
                 T(  31, -31) , T(  32, -32) , T( -29,  29) , T( -30,  30) ,
                 T(  33,  33) , T(  34,  34) , T(  35,  35) , T(  36,  36) ,
                 T(  37,  37) , T(  38,  38) , T(  39,  39) , T(  40,  40) ,
                 T(  41,  41) , T(  42,  42) , T(  43,  43) , T(  44,  44) ,
                 T(  45,  45) , T(  46,  46) , T(  47,  47) , T(  48,  48) ,
                 T(  49,  49) , T(  50,  50) , T(  51,  51) , T(  52,  52) ,
                 T(  53,  53) , T(  54,  54) , T(  55,  55) , T(  56,  56) ,
                 T(  57,  57) , T(  58,  58) , T(  59,  59) , T(  60,  60) ,
                 T(  61,  61) , T(  62,  62) , T(  63,  63) , T(  64,  64) ,
                 T(  67, -67) , T(  68, -68) , T( -65,  65) , T( -66,  66) ,
                 T(  71, -71) , T(  72, -72) , T( -69,  69) , T( -70,  70) ,
                 T(  75, -75) , T(  76, -76) , T( -73,  73) , T( -74,  74) ,
                 T(  79, -79) , T(  80, -80) , T( -77,  77) , T( -78,  78) ,
                 T(  83, -83) , T(  84, -84) , T( -81,  81) , T( -82,  82) ,
                 T(  87, -87) , T(  88, -88) , T( -85,  85) , T( -86,  86) ,
                 T(  91, -91) , T(  92, -92) , T( -89,  89) , T( -90,  90) ,
                 T(  95, -95) , T(  96, -96) , T( -93,  93) , T( -94,  94) ,
                 T(  97,  97) , T(  98,  98) , T(  99,  99) , T( 100, 100) ,
                 T( 101, 101) , T( 102, 102) , T( 103, 103) , T( 104, 104) ,
                 T( 105, 105) , T( 106, 106) , T( 107, 107) , T( 108, 108) ,
                 T( 109, 109) , T( 110, 110) , T( 111, 111) , T( 112, 112) ,
                 T( 113, 113) , T( 114, 114) , T( 115, 115) , T( 116, 116) ,
                 T( 117, 117) , T( 118, 118) , T( 119, 119) , T( 120, 120) ,
                 T( 121, 121) , T( 122, 122) , T( 123, 123) , T( 124, 124) ,
                 T( 125, 125) , T( 126, 126) , T( 127, 127) , T( 128, 128) ,
                 T( 131,-131) , T( 132,-132) , T(-129, 129) , T(-130, 130) ,
                 T( 135,-135) , T( 136,-136) , T(-133, 133) , T(-134, 134) ,
                 T( 139,-139) , T( 140,-140) , T(-137, 137) , T(-138, 138) ,
                 T( 143,-143) , T( 144,-144) , T(-141, 141) , T(-142, 142) ,
                 T( 147,-147) , T( 148,-148) , T(-145, 145) , T(-146, 146) ,
                 T( 151,-151) , T( 152,-152) , T(-149, 149) , T(-150, 150) ,
                 T( 155,-155) , T( 156,-156) , T(-153, 153) , T(-154, 154) ,
                 T( 159,-159) , T( 160,-160) , T(-157, 157) , T(-158, 158) ,
                 T( 161, 161) , T( 162, 162) , T( 163, 163) , T( 164, 164) ,
                 T( 165, 165) , T( 166, 166) , T( 167, 167) , T( 168, 168) ,
                 T( 169, 169) , T( 170, 170) , T( 171, 171) , T( 172, 172) ,
                 T( 173, 173) , T( 174, 174) , T( 175, 175) , T( 176, 176) ,
                 T( 177, 177) , T( 178, 178) , T( 179, 179) , T( 180, 180) ,
                 T( 181, 181) , T( 182, 182) , T( 183, 183) , T( 184, 184) ,
                 T( 185, 185) , T( 186, 186) , T( 187, 187) , T( 188, 188) ,
                 T( 189, 189) , T( 190, 190) , T( 191, 191) , T( 192, 192) ,
                 T( 195,-195) , T( 196,-196) , T(-193, 193) , T(-194, 194) ,
                 T( 199,-199) , T( 200,-200) , T(-197, 197) , T(-198, 198) ,
                 T( 203,-203) , T( 204,-204) , T(-201, 201) , T(-202, 202) ,
                 T( 207,-207) , T( 208,-208) , T(-205, 205) , T(-206, 206) ,
                 T( 211,-211) , T( 212,-212) , T(-209, 209) , T(-210, 210) ,
                 T( 215,-215) , T( 216,-216) , T(-213, 213) , T(-214, 214) ,
                 T( 219,-219) , T( 220,-220) , T(-217, 217) , T(-218, 218) ,
                 T( 223,-223) , T( 224,-224) , T(-221, 221) , T(-222, 222) ,
                 T( 225, 225) , T( 226, 226) , T( 227, 227) , T( 228, 228) ,
                 T( 229, 229) , T( 230, 230) , T( 231, 231) , T( 232, 232) ,
                 T( 233, 233) , T( 234, 234) , T( 235, 235) , T( 236, 236) ,
                 T( 237, 237) , T( 238, 238) , T( 239, 239) , T( 240, 240) ,
                 T( 241, 241) , T( 242, 242) , T( 243, 243) , T( 244, 244) ,
                 T( 245, 245) , T( 246, 246) , T( 247, 247) , T( 248, 248) ,
                 T( 249, 249) , T( 250, 250) , T( 251, 251) , T( 252, 252) ,
                 T( 253, 253) , T( 254, 254) , T( 255, 255) , T( 256, 256) } ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cy26.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cy26.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cy26.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cy26.apply( qclab::Op::Trans , 8 , vec8 ) ;
    check8 = { T(  -3,   3) , T(  -4,   4) , T(   1,  -1) , T(   2,  -2) ,
               T(  -7,   7) , T(  -8,   8) , T(   5,  -5) , T(   6,  -6) ,
               T( -11,  11) , T( -12,  12) , T(   9,  -9) , T(  10, -10) ,
               T( -15,  15) , T( -16,  16) , T(  13, -13) , T(  14, -14) ,
               T( -19,  19) , T( -20,  20) , T(  17, -17) , T(  18, -18) ,
               T( -23,  23) , T( -24,  24) , T(  21, -21) , T(  22, -22) ,
               T( -27,  27) , T( -28,  28) , T(  25, -25) , T(  26, -26) ,
               T( -31,  31) , T( -32,  32) , T(  29, -29) , T(  30, -30) ,
               T(  33,  33) , T(  34,  34) , T(  35,  35) , T(  36,  36) ,
               T(  37,  37) , T(  38,  38) , T(  39,  39) , T(  40,  40) ,
               T(  41,  41) , T(  42,  42) , T(  43,  43) , T(  44,  44) ,
               T(  45,  45) , T(  46,  46) , T(  47,  47) , T(  48,  48) ,
               T(  49,  49) , T(  50,  50) , T(  51,  51) , T(  52,  52) ,
               T(  53,  53) , T(  54,  54) , T(  55,  55) , T(  56,  56) ,
               T(  57,  57) , T(  58,  58) , T(  59,  59) , T(  60,  60) ,
               T(  61,  61) , T(  62,  62) , T(  63,  63) , T(  64,  64) ,
               T( -67,  67) , T( -68,  68) , T(  65, -65) , T(  66, -66) ,
               T( -71,  71) , T( -72,  72) , T(  69, -69) , T(  70, -70) ,
               T( -75,  75) , T( -76,  76) , T(  73, -73) , T(  74, -74) ,
               T( -79,  79) , T( -80,  80) , T(  77, -77) , T(  78, -78) ,
               T( -83,  83) , T( -84,  84) , T(  81, -81) , T(  82, -82) ,
               T( -87,  87) , T( -88,  88) , T(  85, -85) , T(  86, -86) ,
               T( -91,  91) , T( -92,  92) , T(  89, -89) , T(  90, -90) ,
               T( -95,  95) , T( -96,  96) , T(  93, -93) , T(  94, -94) ,
               T(  97,  97) , T(  98,  98) , T(  99,  99) , T( 100, 100) ,
               T( 101, 101) , T( 102, 102) , T( 103, 103) , T( 104, 104) ,
               T( 105, 105) , T( 106, 106) , T( 107, 107) , T( 108, 108) ,
               T( 109, 109) , T( 110, 110) , T( 111, 111) , T( 112, 112) ,
               T( 113, 113) , T( 114, 114) , T( 115, 115) , T( 116, 116) ,
               T( 117, 117) , T( 118, 118) , T( 119, 119) , T( 120, 120) ,
               T( 121, 121) , T( 122, 122) , T( 123, 123) , T( 124, 124) ,
               T( 125, 125) , T( 126, 126) , T( 127, 127) , T( 128, 128) ,
               T(-131, 131) , T(-132, 132) , T( 129,-129) , T( 130,-130) ,
               T(-135, 135) , T(-136, 136) , T( 133,-133) , T( 134,-134) ,
               T(-139, 139) , T(-140, 140) , T( 137,-137) , T( 138,-138) ,
               T(-143, 143) , T(-144, 144) , T( 141,-141) , T( 142,-142) ,
               T(-147, 147) , T(-148, 148) , T( 145,-145) , T( 146,-146) ,
               T(-151, 151) , T(-152, 152) , T( 149,-149) , T( 150,-150) ,
               T(-155, 155) , T(-156, 156) , T( 153,-153) , T( 154,-154) ,
               T(-159, 159) , T(-160, 160) , T( 157,-157) , T( 158,-158) ,
               T( 161, 161) , T( 162, 162) , T( 163, 163) , T( 164, 164) ,
               T( 165, 165) , T( 166, 166) , T( 167, 167) , T( 168, 168) ,
               T( 169, 169) , T( 170, 170) , T( 171, 171) , T( 172, 172) ,
               T( 173, 173) , T( 174, 174) , T( 175, 175) , T( 176, 176) ,
               T( 177, 177) , T( 178, 178) , T( 179, 179) , T( 180, 180) ,
               T( 181, 181) , T( 182, 182) , T( 183, 183) , T( 184, 184) ,
               T( 185, 185) , T( 186, 186) , T( 187, 187) , T( 188, 188) ,
               T( 189, 189) , T( 190, 190) , T( 191, 191) , T( 192, 192) ,
               T(-195, 195) , T(-196, 196) , T( 193,-193) , T( 194,-194) ,
               T(-199, 199) , T(-200, 200) , T( 197,-197) , T( 198,-198) ,
               T(-203, 203) , T(-204, 204) , T( 201,-201) , T( 202,-202) ,
               T(-207, 207) , T(-208, 208) , T( 205,-205) , T( 206,-206) ,
               T(-211, 211) , T(-212, 212) , T( 209,-209) , T( 210,-210) ,
               T(-215, 215) , T(-216, 216) , T( 213,-213) , T( 214,-214) ,
               T(-219, 219) , T(-220, 220) , T( 217,-217) , T( 218,-218) ,
               T(-223, 223) , T(-224, 224) , T( 221,-221) , T( 222,-222) ,
               T( 225, 225) , T( 226, 226) , T( 227, 227) , T( 228, 228) ,
               T( 229, 229) , T( 230, 230) , T( 231, 231) , T( 232, 232) ,
               T( 233, 233) , T( 234, 234) , T( 235, 235) , T( 236, 236) ,
               T( 237, 237) , T( 238, 238) , T( 239, 239) , T( 240, 240) ,
               T( 241, 241) , T( 242, 242) , T( 243, 243) , T( 244, 244) ,
               T( 245, 245) , T( 246, 246) , T( 247, 247) , T( 248, 248) ,
               T( 249, 249) , T( 250, 250) , T( 251, 251) , T( 252, 252) ,
               T( 253, 253) , T( 254, 254) , T( 255, 255) , T( 256, 256) } ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cy26.apply_device( qclab::Op::Trans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

  // control > target
  {
    qclab::qgates::CY< T >  cy10( 1 , 0 , 0 ) ;
    qclab::qgates::CY< T >  cy21( 2 , 1 , 0 ) ;
    qclab::qgates::CY< T >  cy20( 2 , 0 , 0 ) ;
    qclab::qgates::CY< T >  cy31( 3 , 1 , 0 ) ;
    qclab::qgates::CY< T >  cy62( 6 , 2 , 0 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cy10.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { T(3,-3) , T(2,2) , T(-1,1) , T(4,4) } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cy10.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cy10.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cy10.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cy10.apply( qclab::Op::Trans , 2 , vec2 ) ;
    check2 = { T(-3,3) , T(2,2) , T(1,-1) , T(4,4) } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cy10.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    cy10.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { T( 5,-5) , T( 6,-6) , T( 3, 3) , T( 4, 4) ,
                 T(-1, 1) , T(-2, 2) , T( 7, 7) , T( 8, 8) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy10.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy10.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy10.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy10.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { T(-5, 5) , T(-6, 6) , T( 3, 3) , T( 4, 4) ,
               T( 1,-1) , T( 2,-2) , T( 7, 7) , T( 8, 8) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy10.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy21.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T( 3,-3) , T( 2, 2) , T(-1, 1) , T( 4, 4) ,
               T( 7,-7) , T( 6, 6) , T(-5, 5) , T( 8, 8) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy21.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy21.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy21.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy21.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { T(-3, 3) , T( 2, 2) , T( 1,-1) , T( 4, 4) ,
               T(-7, 7) , T( 6, 6) , T( 5,-5) , T( 8, 8) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy21.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy20.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T( 5,-5) , T( 2, 2) , T( 7,-7) , T( 4, 4) ,
               T(-1, 1) , T( 6, 6) , T(-3, 3) , T( 8, 8) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy20.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy20.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy20.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cy20.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { T(-5, 5) , T( 2, 2) , T(-7, 7) , T( 4, 4) ,
               T( 1,-1) , T( 6, 6) , T( 3,-3) , T( 8, 8) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cy20.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    cy21.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { T(  5, -5) , T(  6, -6) , T(  3,  3) , T(  4,  4) ,
                 T( -1,  1) , T( -2,  2) , T(  7,  7) , T(  8,  8) ,
                 T( 13,-13) , T( 14,-14) , T( 11, 11) , T( 12, 12) ,
                 T( -9,  9) , T(-10, 10) , T( 15, 15) , T( 16, 16) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cy21.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cy21.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cy21.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cy21.apply( qclab::Op::Trans , 4 , vec4 ) ;
    check4 = { T( -5,  5) , T( -6,  6) , T(  3,  3) , T(  4,  4) ,
               T(  1, -1) , T(  2, -2) , T(  7,  7) , T(  8,  8) ,
               T(-13, 13) , T(-14, 14) , T( 11, 11) , T( 12, 12) ,
               T(  9, -9) , T( 10,-10) , T( 15, 15) , T( 16, 16) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cy21.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    cy31.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = { T(  9, -9) , T( 10,-10) , T(  3,  3) , T(  4,  4) ,
                 T( 13,-13) , T( 14,-14) , T(  7,  7) , T(  8,  8) ,
                 T( -1,  1) , T( -2,  2) , T( 11, 11) , T( 12, 12) ,
                 T( -5,  5) , T( -6,  6) , T( 15, 15) , T( 16, 16) ,
                 T( 25,-25) , T( 26,-26) , T( 19, 19) , T( 20, 20) ,
                 T( 29,-29) , T( 30,-30) , T( 23, 23) , T( 24, 24) ,
                 T(-17, 17) , T(-18, 18) , T( 27, 27) , T( 28, 28) ,
                 T(-21, 21) , T(-22, 22) , T( 31, 31) , T( 32, 32) } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cy31.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cy31.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cy31.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cy31.apply( qclab::Op::Trans , 5 , vec5 ) ;
    check5 = { T( -9,  9) , T(-10, 10) , T(  3,  3) , T(  4,  4) ,
               T(-13, 13) , T(-14, 14) , T(  7,  7) , T(  8,  8) ,
               T(  1, -1) , T(  2, -2) , T( 11, 11) , T( 12, 12) ,
               T(  5, -5) , T(  6, -6) , T( 15, 15) , T( 16, 16) ,
               T(-25, 25) , T(-26, 26) , T( 19, 19) , T( 20, 20) ,
               T(-29, 29) , T(-30, 30) , T( 23, 23) , T( 24, 24) ,
               T( 17,-17) , T( 18,-18) , T( 27, 27) , T( 28, 28) ,
               T( 21,-21) , T( 22,-22) , T( 31, 31) , T( 32, 32) } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cy31.apply_device( qclab::Op::Trans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    cy62.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    V check8 = { T(  33, -33) , T(  34, -34) , T(   3,   3) , T(   4,   4) ,
                 T(  37, -37) , T(  38, -38) , T(   7,   7) , T(   8,   8) ,
                 T(  41, -41) , T(  42, -42) , T(  11,  11) , T(  12,  12) ,
                 T(  45, -45) , T(  46, -46) , T(  15,  15) , T(  16,  16) ,
                 T(  49, -49) , T(  50, -50) , T(  19,  19) , T(  20,  20) ,
                 T(  53, -53) , T(  54, -54) , T(  23,  23) , T(  24,  24) ,
                 T(  57, -57) , T(  58, -58) , T(  27,  27) , T(  28,  28) ,
                 T(  61, -61) , T(  62, -62) , T(  31,  31) , T(  32,  32) ,
                 T(  -1,   1) , T(  -2,   2) , T(  35,  35) , T(  36,  36) ,
                 T(  -5,   5) , T(  -6,   6) , T(  39,  39) , T(  40,  40) ,
                 T(  -9,   9) , T( -10,  10) , T(  43,  43) , T(  44,  44) ,
                 T( -13,  13) , T( -14,  14) , T(  47,  47) , T(  48,  48) ,
                 T( -17,  17) , T( -18,  18) , T(  51,  51) , T(  52,  52) ,
                 T( -21,  21) , T( -22,  22) , T(  55,  55) , T(  56,  56) ,
                 T( -25,  25) , T( -26,  26) , T(  59,  59) , T(  60,  60) ,
                 T( -29,  29) , T( -30,  30) , T(  63,  63) , T(  64,  64) ,
                 T(  97, -97) , T(  98, -98) , T(  67,  67) , T(  68,  68) ,
                 T( 101,-101) , T( 102,-102) , T(  71,  71) , T(  72,  72) ,
                 T( 105,-105) , T( 106,-106) , T(  75,  75) , T(  76,  76) ,
                 T( 109,-109) , T( 110,-110) , T(  79,  79) , T(  80,  80) ,
                 T( 113,-113) , T( 114,-114) , T(  83,  83) , T(  84,  84) ,
                 T( 117,-117) , T( 118,-118) , T(  87,  87) , T(  88,  88) ,
                 T( 121,-121) , T( 122,-122) , T(  91,  91) , T(  92,  92) ,
                 T( 125,-125) , T( 126,-126) , T(  95,  95) , T(  96,  96) ,
                 T( -65,  65) , T( -66,  66) , T(  99,  99) , T( 100, 100) ,
                 T( -69,  69) , T( -70,  70) , T( 103, 103) , T( 104, 104) ,
                 T( -73,  73) , T( -74,  74) , T( 107, 107) , T( 108, 108) ,
                 T( -77,  77) , T( -78,  78) , T( 111, 111) , T( 112, 112) ,
                 T( -81,  81) , T( -82,  82) , T( 115, 115) , T( 116, 116) ,
                 T( -85,  85) , T( -86,  86) , T( 119, 119) , T( 120, 120) ,
                 T( -89,  89) , T( -90,  90) , T( 123, 123) , T( 124, 124) ,
                 T( -93,  93) , T( -94,  94) , T( 127, 127) , T( 128, 128) ,
                 T( 161,-161) , T( 162,-162) , T( 131, 131) , T( 132, 132) ,
                 T( 165,-165) , T( 166,-166) , T( 135, 135) , T( 136, 136) ,
                 T( 169,-169) , T( 170,-170) , T( 139, 139) , T( 140, 140) ,
                 T( 173,-173) , T( 174,-174) , T( 143, 143) , T( 144, 144) ,
                 T( 177,-177) , T( 178,-178) , T( 147, 147) , T( 148, 148) ,
                 T( 181,-181) , T( 182,-182) , T( 151, 151) , T( 152, 152) ,
                 T( 185,-185) , T( 186,-186) , T( 155, 155) , T( 156, 156) ,
                 T( 189,-189) , T( 190,-190) , T( 159, 159) , T( 160, 160) ,
                 T(-129, 129) , T(-130, 130) , T( 163, 163) , T( 164, 164) ,
                 T(-133, 133) , T(-134, 134) , T( 167, 167) , T( 168, 168) ,
                 T(-137, 137) , T(-138, 138) , T( 171, 171) , T( 172, 172) ,
                 T(-141, 141) , T(-142, 142) , T( 175, 175) , T( 176, 176) ,
                 T(-145, 145) , T(-146, 146) , T( 179, 179) , T( 180, 180) ,
                 T(-149, 149) , T(-150, 150) , T( 183, 183) , T( 184, 184) ,
                 T(-153, 153) , T(-154, 154) , T( 187, 187) , T( 188, 188) ,
                 T(-157, 157) , T(-158, 158) , T( 191, 191) , T( 192, 192) ,
                 T( 225,-225) , T( 226,-226) , T( 195, 195) , T( 196, 196) ,
                 T( 229,-229) , T( 230,-230) , T( 199, 199) , T( 200, 200) ,
                 T( 233,-233) , T( 234,-234) , T( 203, 203) , T( 204, 204) ,
                 T( 237,-237) , T( 238,-238) , T( 207, 207) , T( 208, 208) ,
                 T( 241,-241) , T( 242,-242) , T( 211, 211) , T( 212, 212) ,
                 T( 245,-245) , T( 246,-246) , T( 215, 215) , T( 216, 216) ,
                 T( 249,-249) , T( 250,-250) , T( 219, 219) , T( 220, 220) ,
                 T( 253,-253) , T( 254,-254) , T( 223, 223) , T( 224, 224) ,
                 T(-193, 193) , T(-194, 194) , T( 227, 227) , T( 228, 228) ,
                 T(-197, 197) , T(-198, 198) , T( 231, 231) , T( 232, 232) ,
                 T(-201, 201) , T(-202, 202) , T( 235, 235) , T( 236, 236) ,
                 T(-205, 205) , T(-206, 206) , T( 239, 239) , T( 240, 240) ,
                 T(-209, 209) , T(-210, 210) , T( 243, 243) , T( 244, 244) ,
                 T(-213, 213) , T(-214, 214) , T( 247, 247) , T( 248, 248) ,
                 T(-217, 217) , T(-218, 218) , T( 251, 251) , T( 252, 252) ,
                 T(-221, 221) , T(-222, 222) , T( 255, 255) , T( 256, 256) } ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cy62.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cy62.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cy62.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cy62.apply( qclab::Op::Trans , 8 , vec8 ) ;
    check8 = { T( -33,  33) , T( -34,  34) , T(   3,   3) , T(   4,   4) ,
               T( -37,  37) , T( -38,  38) , T(   7,   7) , T(   8,   8) ,
               T( -41,  41) , T( -42,  42) , T(  11,  11) , T(  12,  12) ,
               T( -45,  45) , T( -46,  46) , T(  15,  15) , T(  16,  16) ,
               T( -49,  49) , T( -50,  50) , T(  19,  19) , T(  20,  20) ,
               T( -53,  53) , T( -54,  54) , T(  23,  23) , T(  24,  24) ,
               T( -57,  57) , T( -58,  58) , T(  27,  27) , T(  28,  28) ,
               T( -61,  61) , T( -62,  62) , T(  31,  31) , T(  32,  32) ,
               T(   1,  -1) , T(   2,  -2) , T(  35,  35) , T(  36,  36) ,
               T(   5,  -5) , T(   6,  -6) , T(  39,  39) , T(  40,  40) ,
               T(   9,  -9) , T(  10, -10) , T(  43,  43) , T(  44,  44) ,
               T(  13, -13) , T(  14, -14) , T(  47,  47) , T(  48,  48) ,
               T(  17, -17) , T(  18, -18) , T(  51,  51) , T(  52,  52) ,
               T(  21, -21) , T(  22, -22) , T(  55,  55) , T(  56,  56) ,
               T(  25, -25) , T(  26, -26) , T(  59,  59) , T(  60,  60) ,
               T(  29, -29) , T(  30, -30) , T(  63,  63) , T(  64,  64) ,
               T( -97,  97) , T( -98,  98) , T(  67,  67) , T(  68,  68) ,
               T(-101, 101) , T(-102, 102) , T(  71,  71) , T(  72,  72) ,
               T(-105, 105) , T(-106, 106) , T(  75,  75) , T(  76,  76) ,
               T(-109, 109) , T(-110, 110) , T(  79,  79) , T(  80,  80) ,
               T(-113, 113) , T(-114, 114) , T(  83,  83) , T(  84,  84) ,
               T(-117, 117) , T(-118, 118) , T(  87,  87) , T(  88,  88) ,
               T(-121, 121) , T(-122, 122) , T(  91,  91) , T(  92,  92) ,
               T(-125, 125) , T(-126, 126) , T(  95,  95) , T(  96,  96) ,
               T(  65, -65) , T(  66, -66) , T(  99,  99) , T( 100, 100) ,
               T(  69, -69) , T(  70, -70) , T( 103, 103) , T( 104, 104) ,
               T(  73, -73) , T(  74, -74) , T( 107, 107) , T( 108, 108) ,
               T(  77, -77) , T(  78, -78) , T( 111, 111) , T( 112, 112) ,
               T(  81, -81) , T(  82, -82) , T( 115, 115) , T( 116, 116) ,
               T(  85, -85) , T(  86, -86) , T( 119, 119) , T( 120, 120) ,
               T(  89, -89) , T(  90, -90) , T( 123, 123) , T( 124, 124) ,
               T(  93, -93) , T(  94, -94) , T( 127, 127) , T( 128, 128) ,
               T(-161, 161) , T(-162, 162) , T( 131, 131) , T( 132, 132) ,
               T(-165, 165) , T(-166, 166) , T( 135, 135) , T( 136, 136) ,
               T(-169, 169) , T(-170, 170) , T( 139, 139) , T( 140, 140) ,
               T(-173, 173) , T(-174, 174) , T( 143, 143) , T( 144, 144) ,
               T(-177, 177) , T(-178, 178) , T( 147, 147) , T( 148, 148) ,
               T(-181, 181) , T(-182, 182) , T( 151, 151) , T( 152, 152) ,
               T(-185, 185) , T(-186, 186) , T( 155, 155) , T( 156, 156) ,
               T(-189, 189) , T(-190, 190) , T( 159, 159) , T( 160, 160) ,
               T( 129,-129) , T( 130,-130) , T( 163, 163) , T( 164, 164) ,
               T( 133,-133) , T( 134,-134) , T( 167, 167) , T( 168, 168) ,
               T( 137,-137) , T( 138,-138) , T( 171, 171) , T( 172, 172) ,
               T( 141,-141) , T( 142,-142) , T( 175, 175) , T( 176, 176) ,
               T( 145,-145) , T( 146,-146) , T( 179, 179) , T( 180, 180) ,
               T( 149,-149) , T( 150,-150) , T( 183, 183) , T( 184, 184) ,
               T( 153,-153) , T( 154,-154) , T( 187, 187) , T( 188, 188) ,
               T( 157,-157) , T( 158,-158) , T( 191, 191) , T( 192, 192) ,
               T(-225, 225) , T(-226, 226) , T( 195, 195) , T( 196, 196) ,
               T(-229, 229) , T(-230, 230) , T( 199, 199) , T( 200, 200) ,
               T(-233, 233) , T(-234, 234) , T( 203, 203) , T( 204, 204) ,
               T(-237, 237) , T(-238, 238) , T( 207, 207) , T( 208, 208) ,
               T(-241, 241) , T(-242, 242) , T( 211, 211) , T( 212, 212) ,
               T(-245, 245) , T(-246, 246) , T( 215, 215) , T( 216, 216) ,
               T(-249, 249) , T(-250, 250) , T( 219, 219) , T( 220, 220) ,
               T(-253, 253) , T(-254, 254) , T( 223, 223) , T( 224, 224) ,
               T( 193,-193) , T( 194,-194) , T( 227, 227) , T( 228, 228) ,
               T( 197,-197) , T( 198,-198) , T( 231, 231) , T( 232, 232) ,
               T( 201,-201) , T( 202,-202) , T( 235, 235) , T( 236, 236) ,
               T( 205,-205) , T( 206,-206) , T( 239, 239) , T( 240, 240) ,
               T( 209,-209) , T( 210,-210) , T( 243, 243) , T( 244, 244) ,
               T( 213,-213) , T( 214,-214) , T( 247, 247) , T( 248, 248) ,
               T( 217,-217) , T( 218,-218) , T( 251, 251) , T( 252, 252) ,
               T( 221,-221) , T( 222,-222) , T( 255, 255) , T( 256, 256) } ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cy62.apply_device( qclab::Op::Trans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

}


/*
 * complex float
 */
TEST( qclab_qgates_CY , complex_float ) {
  test_qclab_qgates_CY< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_CY , complex_double ) {
  test_qclab_qgates_CY< std::complex< double > >() ;
}

