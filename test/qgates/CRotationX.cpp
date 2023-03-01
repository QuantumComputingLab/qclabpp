#include <gtest/gtest.h>
#include "qclab/qgates/CRotationX.hpp"
#include <numeric>

template <typename T>
void test_qclab_qgates_CRotationX() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R tol = 10 * std::numeric_limits< R >::epsilon() ;

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

    // rotation, theta, sin, cos
    qclab::QRotation< R > rot ;
    EXPECT_TRUE( crotx.rotation() == rot ) ;
    EXPECT_EQ( crotx.theta() , 0 ) ;
    EXPECT_EQ( crotx.cos() , 1 ) ;
    EXPECT_EQ( crotx.sin() , 0 ) ;

    // update(rot)
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    rot.update( pi/2 ) ;
    crotx.update( rot ) ;
    EXPECT_NEAR( crotx.theta() , pi/2 , tol ) ;
    EXPECT_NEAR( crotx.cos() , cos , tol ) ;
    EXPECT_NEAR( crotx.sin() , sin , tol ) ;

    // update(theta)
    crotx.update( pi ) ;
    EXPECT_NEAR( crotx.theta() , pi , tol ) ;
    EXPECT_NEAR( crotx.cos() , 0 , tol ) ;
    EXPECT_NEAR( crotx.sin() , 1 , tol ) ;

    // update(cos,sin)
    crotx.update( cos , sin ) ;
    EXPECT_NEAR( crotx.theta() , pi/2 , tol ) ;
    EXPECT_NEAR( crotx.cos() , cos , tol ) ;
    EXPECT_NEAR( crotx.sin() , sin , tol ) ;
  }

  {
    const R theta = pi/4 ;
    qclab::QRotation< R >  rot( theta ) ;
    qclab::qgates::CRotationX< T >  crotx( 1 , 3 , rot ) ;

    EXPECT_EQ( crotx.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( crotx.fixed() ) ;          // fixed
    EXPECT_TRUE( crotx.controlled() ) ;      // controlled
    EXPECT_EQ( crotx.control() , 1 ) ;       // control
    EXPECT_EQ( crotx.target() , 3 ) ;        // target
    EXPECT_EQ( crotx.controlState() , 1 ) ;  // controlState

    EXPECT_NEAR( crotx.theta() , theta , tol ) ;
    EXPECT_NEAR( crotx.cos() , std::cos( theta/2 ) , tol ) ;
    EXPECT_NEAR( crotx.sin() , std::sin( theta/2 ) , tol ) ;
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

    EXPECT_NEAR( crotx.theta() , theta , tol ) ;
    EXPECT_NEAR( crotx.cos() , std::cos( theta/2 ) , tol ) ;
    EXPECT_NEAR( crotx.sin() , std::sin( theta/2 ) , tol ) ;
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

    EXPECT_NEAR( crotx.theta() , 2*theta , tol ) ;
    EXPECT_NEAR( crotx.cos() , cos , tol ) ;
    EXPECT_NEAR( crotx.sin() , sin , tol ) ;
  }

  using V = std::vector< T > ;

  const T c = T( std::cos( pi/8 ) ) ;
  const T s = T( 0 , std::sin( pi/8 ) ) ;

  const V v2 = { 3 , 5 , 2 , 7 } ;
  const V v3 = { 3 , 5 , 2 , 7 , 4 , 1 , 8 , 3 } ;
  const V v4 = { 3 , 5 , 2 , 7 , 4 , 1 , 8 , 3 , 7 , 2 , 5 , 6 , 8 , 9 , 5 , 1};
  V vt( 32 ) ; std::iota( vt.begin() , vt.end() , 0 ) ;
  const V v5 = vt ;
  vt = V( 256 ) ; std::iota( vt.begin() , vt.end() , 0 ) ;
  const V v8 = vt ;

  //
  // CNOTs |1> controlled
  //

  // control < target
  {
    qclab::qgates::CRotationX< T >  crotx01( 0 , 1 , pi/4 ) ;
    qclab::qgates::CRotationX< T >  crotx12( 1 , 2 , pi/4 ) ;
    qclab::qgates::CRotationX< T >  crotx02( 0 , 2 , pi/4 ) ;
    qclab::qgates::CRotationX< T >  crotx13( 1 , 3 , pi/4 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    crotx01.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 3 , 5 , T(2)*c - T(7)*s , -T(2)*s + T(7)*c } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { crotx01.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    crotx01.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check2 = { 3 , 5 , T(2)*c + T(7)*s , T(2)*s + T(7)*c } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { crotx01.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    crotx01.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 3 , 5 , 2 , 7 , T(4)*c - T(8)*s , T(1)*c - T(3)*s ,
                                -T(4)*s + T(8)*c ,-T(1)*s + T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx01.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx01.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , 7 , T(4)*c + T(8)*s , T(1)*c + T(3)*s ,
                               T(4)*s + T(8)*c , T(1)*s + T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx01.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx12.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , T(2)*c - T(7)*s , -T(2)*s + T(7)*c ,
               4 , 1 , T(8)*c - T(3)*s , -T(8)*s + T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx12.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx12.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , T(2)*c + T(7)*s , T(2)*s + T(7)*c ,
               4 , 1 , T(8)*c + T(3)*s , T(8)*s + T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx12.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx02.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , 7 , T(4)*c - T(1)*s , -T(4)*s + T(1)*c ,
                               T(8)*c - T(3)*s , -T(8)*s + T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx02.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx02.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , 7 , T(4)*c + T(1)*s , T(4)*s + T(1)*c ,
                               T(8)*c + T(3)*s , T(8)*s + T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx02.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    crotx12.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { 3 , 5 , 2 , 7 , T(4)*c - T(8)*s , T(1)*c - T(3)*s ,
                                -T(4)*s + T(8)*c ,-T(1)*s + T(3)*c ,
                 7 , 2 , 5 , 6 , T(8)*c - T(5)*s , T(9)*c - T(1)*s ,
                                -T(8)*s + T(5)*c ,-T(9)*s + T(1)*c } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { crotx12.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    crotx12.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { 3 , 5 , 2 , 7 , T(4)*c + T(8)*s , T(1)*c + T(3)*s ,
                               T(4)*s + T(8)*c , T(1)*s + T(3)*c ,
               7 , 2 , 5 , 6 , T(8)*c + T(5)*s , T(9)*c + T(1)*s ,
                               T(8)*s + T(5)*c , T(9)*s + T(1)*c } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { crotx12.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    crotx13.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = {  0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,
                 T( 8)*c - T(10)*s , T( 9)*c - T(11)*s ,
                -T( 8)*s + T(10)*c ,-T( 9)*s + T(11)*c ,
                 T(12)*c - T(14)*s , T(13)*c - T(15)*s ,
                -T(12)*s + T(14)*c ,-T(13)*s + T(15)*c ,
                 16 , 17 , 18 , 19 , 20 , 21 , 22 , 23 ,
                 T(24)*c - T(26)*s , T(25)*c - T(27)*s ,
                -T(24)*s + T(26)*c ,-T(25)*s + T(27)*c ,
                 T(28)*c - T(30)*s , T(29)*c - T(31)*s ,
                -T(28)*s + T(30)*c ,-T(29)*s + T(31)*c } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { crotx13.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    crotx13.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    check5 = {  0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,
               T( 8)*c + T(10)*s , T( 9)*c + T(11)*s ,
               T( 8)*s + T(10)*c , T( 9)*s + T(11)*c ,
               T(12)*c + T(14)*s , T(13)*c + T(15)*s ,
               T(12)*s + T(14)*c , T(13)*s + T(15)*c ,
               16 , 17 , 18 , 19 , 20 , 21 , 22 , 23 ,
               T(24)*c + T(26)*s , T(25)*c + T(27)*s ,
               T(24)*s + T(26)*c , T(25)*s + T(27)*c ,
               T(28)*c + T(30)*s , T(29)*c + T(31)*s ,
               T(28)*s + T(30)*c , T(29)*s + T(31)*c } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { crotx13.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif
  }

  // control > target
  {
    qclab::qgates::CRotationX< T >  crotx10( 1 , 0 , pi/4 ) ;
    qclab::qgates::CRotationX< T >  crotx21( 2 , 1 , pi/4 ) ;
    qclab::qgates::CRotationX< T >  crotx20( 2 , 0 , pi/4 ) ;
    qclab::qgates::CRotationX< T >  crotx31( 3 , 1 , pi/4 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    crotx10.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 3 , T(5)*c - T(7)*s , 2 , -T(5)*s + T(7)*c } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { crotx10.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    crotx10.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check2 = { 3 , T(5)*c + T(7)*s , 2 , T(5)*s + T(7)*c } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { crotx10.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    crotx10.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 3 , 5 , T(2)*c - T(8)*s , T(7)*c - T(3)*s ,
                 4 , 1 ,-T(2)*s + T(8)*c ,-T(7)*s + T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx10.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx10.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , T(2)*c + T(8)*s , T(7)*c + T(3)*s ,
               4 , 1 , T(2)*s + T(8)*c , T(7)*s + T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx10.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx21.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , T(5)*c - T(7)*s , 2 , -T(5)*s + T(7)*c ,
               4 , T(1)*c - T(3)*s , 8 , -T(1)*s + T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx21.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx21.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , T(5)*c + T(7)*s , 2 , T(5)*s + T(7)*c ,
               4 , T(1)*c + T(3)*s , 8 , T(1)*s + T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx21.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx20.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , T(5)*c - T(1)*s , 2 , T(7)*c - T(3)*s ,
               4 ,-T(5)*s + T(1)*c , 8 ,-T(7)*s + T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx20.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx20.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , T(5)*c + T(1)*s , 2 , T(7)*c + T(3)*s ,
               4 , T(5)*s + T(1)*c , 8 , T(7)*s + T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx20.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    crotx21.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { 3 , 5 , T(2)*c - T(8)*s , T(7)*c - T(3)*s ,
                 4 , 1 ,-T(2)*s + T(8)*c ,-T(7)*s + T(3)*c ,
                 7 , 2 , T(5)*c - T(5)*s , T(6)*c - T(1)*s ,
                 8 , 9 ,-T(5)*s + T(5)*c ,-T(6)*s + T(1)*c } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { crotx21.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    crotx21.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { 3 , 5 , T(2)*c + T(8)*s , T(7)*c + T(3)*s ,
               4 , 1 , T(2)*s + T(8)*c , T(7)*s + T(3)*c ,
               7 , 2 , T(5)*c + T(5)*s , T(6)*c + T(1)*s ,
               8 , 9 , T(5)*s + T(5)*c , T(6)*s + T(1)*c } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { crotx21.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    crotx31.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = {  0 ,  1 , T( 2)*c - T(10)*s , T( 3)*c - T(11)*s ,
                  4 ,  5 , T( 6)*c - T(14)*s , T( 7)*c - T(15)*s ,
                  8 ,  9 ,-T( 2)*s + T(10)*c ,-T( 3)*s + T(11)*c ,
                 12 , 13 ,-T( 6)*s + T(14)*c ,-T( 7)*s + T(15)*c ,
                 16 , 17 , T(18)*c - T(26)*s , T(19)*c - T(27)*s ,
                 20 , 21 , T(22)*c - T(30)*s , T(23)*c - T(31)*s ,
                 24 , 25 ,-T(18)*s + T(26)*c ,-T(19)*s + T(27)*c ,
                 28 , 29 ,-T(22)*s + T(30)*c ,-T(23)*s + T(31)*c } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { crotx31.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    crotx31.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    check5 = {  0 ,  1 , T( 2)*c + T(10)*s , T( 3)*c + T(11)*s ,
                4 ,  5 , T( 6)*c + T(14)*s , T( 7)*c + T(15)*s ,
                8 ,  9 , T( 2)*s + T(10)*c , T( 3)*s + T(11)*c ,
               12 , 13 , T( 6)*s + T(14)*c , T( 7)*s + T(15)*c ,
               16 , 17 , T(18)*c + T(26)*s , T(19)*c + T(27)*s ,
               20 , 21 , T(22)*c + T(30)*s , T(23)*c + T(31)*s ,
               24 , 25 , T(18)*s + T(26)*c , T(19)*s + T(27)*c ,
               28 , 29 , T(22)*s + T(30)*c , T(23)*s + T(31)*c } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { crotx31.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif
  }

  //
  // CNOTs |0> controlled
  //

  // control < target
  {
    qclab::qgates::CRotationX< T >  crotx01( 0 , 1 , pi/4 , 0 ) ;
    qclab::qgates::CRotationX< T >  crotx12( 1 , 2 , pi/4 , 0 ) ;
    qclab::qgates::CRotationX< T >  crotx02( 0 , 2 , pi/4 , 0 ) ;
    qclab::qgates::CRotationX< T >  crotx13( 1 , 3 , pi/4 , 0 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    crotx01.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { T(3)*c - T(5)*s , -T(3)*s + T(5)*c , 2 , 7 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { crotx01.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    crotx01.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check2 = { T(3)*c + T(5)*s , T(3)*s + T(5)*c , 2 , 7 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { crotx01.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    crotx01.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { T(3)*c - T(2)*s , T(5)*c - T(7)*s ,
                -T(3)*s + T(2)*c ,-T(5)*s + T(7)*c , 4 , 1 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx01.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx01.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { T(3)*c + T(2)*s , T(5)*c + T(7)*s ,
               T(3)*s + T(2)*c , T(5)*s + T(7)*c , 4 , 1 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx01.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx12.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T(3)*c - T(5)*s , -T(3)*s + T(5)*c , 2 , 7 ,
               T(4)*c - T(1)*s , -T(4)*s + T(1)*c , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx12.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx12.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { T(3)*c + T(5)*s , T(3)*s + T(5)*c , 2 , 7 ,
               T(4)*c + T(1)*s , T(4)*s + T(1)*c , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx12.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx02.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T(3)*c - T(5)*s , -T(3)*s + T(5)*c ,
               T(2)*c - T(7)*s , -T(2)*s + T(7)*c , 4 , 1 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx02.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx02.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { T(3)*c + T(5)*s , T(3)*s + T(5)*c ,
               T(2)*c + T(7)*s , T(2)*s + T(7)*c , 4 , 1 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx02.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    crotx12.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { T(3)*c - T(2)*s , T(5)*c - T(7)*s ,
                -T(3)*s + T(2)*c ,-T(5)*s + T(7)*c , 4 , 1 , 8 , 3 ,
                 T(7)*c - T(5)*s , T(2)*c - T(6)*s ,
                -T(7)*s + T(5)*c ,-T(2)*s + T(6)*c , 8 , 9 , 5 , 1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { crotx12.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    crotx12.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { T(3)*c + T(2)*s , T(5)*c + T(7)*s ,
               T(3)*s + T(2)*c , T(5)*s + T(7)*c , 4 , 1 , 8 , 3 ,
               T(7)*c + T(5)*s , T(2)*c + T(6)*s ,
               T(7)*s + T(5)*c , T(2)*s + T(6)*c , 8 , 9 , 5 , 1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { crotx12.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    crotx13.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = { T( 0)*c - T( 2)*s , T( 1)*c - T( 3)*s ,
                -T( 0)*s + T( 2)*c ,-T( 1)*s + T( 3)*c ,
                 T( 4)*c - T( 6)*s , T( 5)*c - T( 7)*s ,
                -T( 4)*s + T( 6)*c ,-T( 5)*s + T( 7)*c ,
                  8 ,  9 , 10 , 11 , 12 , 13 , 14 , 15 ,
                 T(16)*c - T(18)*s , T(17)*c - T(19)*s ,
                -T(16)*s + T(18)*c ,-T(17)*s + T(19)*c ,
                 T(20)*c - T(22)*s , T(21)*c - T(23)*s ,
                -T(20)*s + T(22)*c ,-T(21)*s + T(23)*c ,
                 24 , 25 , 26 , 27 , 28 , 29 , 30 , 31 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { crotx13.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    crotx13.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    check5 = { T( 0)*c + T( 2)*s , T( 1)*c + T( 3)*s ,
               T( 0)*s + T( 2)*c , T( 1)*s + T( 3)*c ,
               T( 4)*c + T( 6)*s , T( 5)*c + T( 7)*s ,
               T( 4)*s + T( 6)*c , T( 5)*s + T( 7)*c ,
                8 ,  9 , 10 , 11 , 12 , 13 , 14 , 15 ,
               T(16)*c + T(18)*s , T(17)*c + T(19)*s ,
               T(16)*s + T(18)*c , T(17)*s + T(19)*c ,
               T(20)*c + T(22)*s , T(21)*c + T(23)*s ,
               T(20)*s + T(22)*c , T(21)*s + T(23)*c ,
               24 , 25 , 26 , 27 , 28 , 29 , 30 , 31 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { crotx13.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif
  }

  // control > target
  {
    qclab::qgates::CRotationX< T >  crotx10( 1 , 0 , pi/4 , 0 ) ;
    qclab::qgates::CRotationX< T >  crotx21( 2 , 1 , pi/4 , 0 ) ;
    qclab::qgates::CRotationX< T >  crotx20( 2 , 0 , pi/4 , 0 ) ;
    qclab::qgates::CRotationX< T >  crotx31( 3 , 1 , pi/4 , 0 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    crotx10.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { T(3)*c - T(2)*s , 5 , -T(3)*s + T(2)*c , 7 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { crotx10.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    crotx10.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check2 = { T(3)*c + T(2)*s , 5 , T(3)*s + T(2)*c , 7 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { crotx10.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    crotx10.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { T(3)*c - T(4)*s , T(5)*c - T(1)*s , 2 , 7 ,
                -T(3)*s + T(4)*c ,-T(5)*s + T(1)*c , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx10.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx10.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { T(3)*c + T(4)*s , T(5)*c + T(1)*s , 2 , 7 ,
               T(3)*s + T(4)*c , T(5)*s + T(1)*c , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx10.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx21.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T(3)*c - T(2)*s , 5 , -T(3)*s + T(2)*c , 7 ,
               T(4)*c - T(8)*s , 1 , -T(4)*s + T(8)*c , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx21.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx21.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { T(3)*c + T(2)*s , 5 , T(3)*s + T(2)*c , 7 ,
               T(4)*c + T(8)*s , 1 , T(4)*s + T(8)*c , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx21.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx20.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T(3)*c - T(4)*s , 5 , T(2)*c - T(8)*s , 7 ,
              -T(3)*s + T(4)*c , 1 ,-T(2)*s + T(8)*c , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx20.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotx20.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { T(3)*c + T(4)*s , 5 , T(2)*c + T(8)*s , 7 ,
               T(3)*s + T(4)*c , 1 , T(2)*s + T(8)*c , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotx20.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    crotx21.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { T(3)*c - T(4)*s , T(5)*c - T(1)*s , 2 , 7 ,
                -T(3)*s + T(4)*c ,-T(5)*s + T(1)*c , 8 , 3 ,
                 T(7)*c - T(8)*s , T(2)*c - T(9)*s , 5 , 6 ,
                -T(7)*s + T(8)*c ,-T(2)*s + T(9)*c , 5 , 1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { crotx21.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    crotx21.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { T(3)*c + T(4)*s , T(5)*c + T(1)*s , 2 , 7 ,
               T(3)*s + T(4)*c , T(5)*s + T(1)*c , 8 , 3 ,
               T(7)*c + T(8)*s , T(2)*c + T(9)*s , 5 , 6 ,
               T(7)*s + T(8)*c , T(2)*s + T(9)*c , 5 , 1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { crotx21.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    crotx31.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = { T( 0)*c - T( 8)*s , T( 1)*c - T( 9)*s ,  2 ,  3 ,
                 T( 4)*c - T(12)*s , T( 5)*c - T(13)*s ,  6 ,  7 ,
                -T( 0)*s + T( 8)*c ,-T( 1)*s + T( 9)*c , 10 , 11 ,
                -T( 4)*s + T(12)*c ,-T( 5)*s + T(13)*c , 14 , 15 ,
                 T(16)*c - T(24)*s , T(17)*c - T(25)*s , 18 , 19 ,
                 T(20)*c - T(28)*s , T(21)*c - T(29)*s , 22 , 23 ,
                -T(16)*s + T(24)*c ,-T(17)*s + T(25)*c , 26 , 27 ,
                -T(20)*s + T(28)*c ,-T(21)*s + T(29)*c , 30 , 31 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { crotx31.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    crotx31.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    check5 = { T( 0)*c + T( 8)*s , T( 1)*c + T( 9)*s ,  2 ,  3 ,
               T( 4)*c + T(12)*s , T( 5)*c + T(13)*s ,  6 ,  7 ,
               T( 0)*s + T( 8)*c , T( 1)*s + T( 9)*c , 10 , 11 ,
               T( 4)*s + T(12)*c , T( 5)*s + T(13)*c , 14 , 15 ,
               T(16)*c + T(24)*s , T(17)*c + T(25)*s , 18 , 19 ,
               T(20)*c + T(28)*s , T(21)*c + T(29)*s , 22 , 23 ,
               T(16)*s + T(24)*c , T(17)*s + T(25)*c , 26 , 27 ,
               T(20)*s + T(28)*c , T(21)*s + T(29)*c , 30 , 31 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { crotx31.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif
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

