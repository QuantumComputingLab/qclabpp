#include <gtest/gtest.h>
#include "qclab/qgates/CRotationZ.hpp"
#include <numeric>

template <typename T>
void test_qclab_qgates_CRotationZ() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R tol = 10 * std::numeric_limits< R >::epsilon() ;

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
    EXPECT_NEAR( crotz.theta() , pi/2 , tol ) ;
    EXPECT_NEAR( crotz.cos() , cos , tol ) ;
    EXPECT_NEAR( crotz.sin() , sin , tol ) ;

    // update(theta)
    crotz.update( pi ) ;
    EXPECT_NEAR( crotz.theta() , pi , tol ) ;
    EXPECT_NEAR( crotz.cos() , 0 , tol ) ;
    EXPECT_NEAR( crotz.sin() , 1 , tol ) ;

    // update(cos,sin)
    crotz.update( cos , sin ) ;
    EXPECT_NEAR( crotz.theta() , pi/2 , tol ) ;
    EXPECT_NEAR( crotz.cos() , cos , tol ) ;
    EXPECT_NEAR( crotz.sin() , sin , tol ) ;
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

    EXPECT_NEAR( crotz.theta() , theta , tol ) ;
    EXPECT_NEAR( crotz.cos() , std::cos( theta/2 ) , tol ) ;
    EXPECT_NEAR( crotz.sin() , std::sin( theta/2 ) , tol ) ;
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

    EXPECT_NEAR( crotz.theta() , theta , tol ) ;
    EXPECT_NEAR( crotz.cos() , std::cos( theta/2 ) , tol ) ;
    EXPECT_NEAR( crotz.sin() , std::sin( theta/2 ) , tol ) ;
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

    EXPECT_NEAR( crotz.theta() , 2*theta , tol ) ;
    EXPECT_NEAR( crotz.cos() , cos , tol ) ;
    EXPECT_NEAR( crotz.sin() , sin , tol ) ;
  }

  using V = std::vector< T > ;

  const T c1 = T( std::cos( pi/4 ) , -std::sin( pi/4 ) ) ;
  const T c2 = T( std::cos( pi/4 ) ,  std::sin( pi/4 ) ) ;

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
    qclab::qgates::CRotationZ< T >  crotz01( 0 , 1 , pi/2 ) ;
    qclab::qgates::CRotationZ< T >  crotz12( 1 , 2 , pi/2 ) ;
    qclab::qgates::CRotationZ< T >  crotz02( 0 , 2 , pi/2 ) ;
    qclab::qgates::CRotationZ< T >  crotz13( 1 , 3 , pi/2 ) ;
    qclab::qgates::CRotationZ< T >  crotz26( 2 , 6 , pi/2 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    crotz01.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 3 , 5 , T(2)*c1 , T(7)*c2 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { crotz01.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    crotz01.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check2 = { 3 , 5 , T(2)*c2 , T(7)*c1 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { crotz01.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    crotz01.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 3 , 5 , 2 , 7 , T(4)*c1 , T(1)*c1 , T(8)*c2 , T(3)*c2 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz01.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz01.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , 7 , T(4)*c2 , T(1)*c2 , T(8)*c1 , T(3)*c1 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz01.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz12.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , T(2)*c1 , T(7)*c2 , 4 , 1 , T(8)*c1 , T(3)*c2 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz12.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz12.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , T(2)*c2 , T(7)*c1 , 4 , 1 , T(8)*c2 , T(3)*c1 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz12.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz02.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , 7 , T(4)*c1 , T(1)*c2 , T(8)*c1 , T(3)*c2 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz02.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz02.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , 7 , T(4)*c2 , T(1)*c1 , T(8)*c2 , T(3)*c1 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz02.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    crotz12.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { 3 , 5 , 2 , 7 , T(4)*c1 , T(1)*c1 , T(8)*c2 , T(3)*c2 ,
                 7 , 2 , 5 , 6 , T(8)*c1 , T(9)*c1 , T(5)*c2 , T(1)*c2 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { crotz12.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    crotz12.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { 3 , 5 , 2 , 7 , T(4)*c2 , T(1)*c2 , T(8)*c1 , T(3)*c1 ,
               7 , 2 , 5 , 6 , T(8)*c2 , T(9)*c2 , T(5)*c1 , T(1)*c1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { crotz12.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    crotz13.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = {  0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,
                 T( 8)*c1 , T( 9)*c1 , T(10)*c2 , T(11)*c2 ,
                 T(12)*c1 , T(13)*c1 , T(14)*c2 , T(15)*c2 ,
                 16 , 17 , 18 , 19 , 20 , 21 , 22 , 23 ,
                 T(24)*c1 , T(25)*c1 , T(26)*c2 , T(27)*c2 ,
                 T(28)*c1 , T(29)*c1 , T(30)*c2 , T(31)*c2 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { crotz13.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    crotz13.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    check5 = {  0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,
               T( 8)*c2 , T( 9)*c2 , T(10)*c1 , T(11)*c1 ,
               T(12)*c2 , T(13)*c2 , T(14)*c1 , T(15)*c1 ,
               16 , 17 , 18 , 19 , 20 , 21 , 22 , 23 ,
               T(24)*c2 , T(25)*c2 , T(26)*c1 , T(27)*c1 ,
               T(28)*c2 , T(29)*c2 , T(30)*c1 , T(31)*c1 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { crotz13.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    crotz26.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    std::vector< int > i8 = {  34 ,  35 ,  38 ,  39 ,  42 ,  43 ,  46 ,  47 ,
                               50 ,  51 ,  54 ,  55 ,  58 ,  59 ,  62 ,  63 ,
                               98 ,  99 , 102 , 103 , 106 , 107 , 110 , 111 ,
                              114 , 115 , 118 , 119 , 122 , 123 , 126 , 127 ,
                              162 , 163 , 166 , 167 , 170 , 171 , 174 , 175 ,
                              178 , 179 , 182 , 183 , 186 , 187 , 190 , 191 ,
                              226 , 227 , 230 , 231 , 234 , 235 , 238 , 239 ,
                              242 , 243 , 246 , 247 , 250 , 251 , 254 , 255 } ;
    V check8 = v8 ; for ( auto& i: i8 ) { check8[i-2] *= c1 ; check8[i] *= c2 ;}
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { crotz26.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    crotz26.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    check8 = v8 ; for ( auto& i: i8 ) { check8[i-2] *= c2 ; check8[i] *= c1 ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { crotz26.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

  // control > target
  {
    qclab::qgates::CRotationZ< T >  crotz10( 1 , 0 , pi/2 ) ;
    qclab::qgates::CRotationZ< T >  crotz21( 2 , 1 , pi/2 ) ;
    qclab::qgates::CRotationZ< T >  crotz20( 2 , 0 , pi/2 ) ;
    qclab::qgates::CRotationZ< T >  crotz31( 3 , 1 , pi/2 ) ;
    qclab::qgates::CRotationZ< T >  crotz62( 6 , 2 , pi/2 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    crotz10.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 3 , T(5)*c1 , 2 , T(7)*c2 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { crotz10.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    crotz10.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check2 = { 3 , T(5)*c2 , 2 , T(7)*c1 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { crotz10.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    crotz10.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 3 , 5 , T(2)*c1 , T(7)*c1 , 4 , 1 , T(8)*c2 , T(3)*c2 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz10.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz10.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , T(2)*c2 , T(7)*c2 , 4 , 1 , T(8)*c1 , T(3)*c1 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz10.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz21.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , T(5)*c1 , 2 , T(7)*c2 , 4 , T(1)*c1 , 8 , T(3)*c2 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz21.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz21.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , T(5)*c2 , 2 , T(7)*c1 , 4 , T(1)*c2 , 8 , T(3)*c1 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz21.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz20.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , T(5)*c1 , 2 , T(7)*c1 , 4 , T(1)*c2 , 8 , T(3)*c2 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz20.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz20.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , T(5)*c2 , 2 , T(7)*c2 , 4 , T(1)*c1 , 8 , T(3)*c1 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz20.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    crotz21.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { 3 , 5 , T(2)*c1 , T(7)*c1 , 4 , 1 , T(8)*c2 , T(3)*c2 ,
                 7 , 2 , T(5)*c1 , T(6)*c1 , 8 , 9 , T(5)*c2 , T(1)*c2 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { crotz21.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    crotz21.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { 3 , 5 , T(2)*c2 , T(7)*c2 , 4 , 1 , T(8)*c1 , T(3)*c1 ,
               7 , 2 , T(5)*c2 , T(6)*c2 , 8 , 9 , T(5)*c1 , T(1)*c1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { crotz21.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    crotz31.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = {  0 ,  1 , T( 2)*c1 , T( 3)*c1 ,  4 ,  5 , T( 6)*c1 , T( 7)*c1 ,
                  8 ,  9 , T(10)*c2 , T(11)*c2 , 12 , 13 , T(14)*c2 , T(15)*c2 ,
                 16 , 17 , T(18)*c1 , T(19)*c1 , 20 , 21 , T(22)*c1 , T(23)*c1 ,
                 24 , 25 , T(26)*c2 , T(27)*c2 , 28 , 29 , T(30)*c2 , T(31)*c2};
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { crotz31.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    crotz31.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    check5 = {  0 ,  1 , T( 2)*c2 , T( 3)*c2 ,  4 ,  5 , T( 6)*c2 , T( 7)*c2 ,
                8 ,  9 , T(10)*c1 , T(11)*c1 , 12 , 13 , T(14)*c1 , T(15)*c1 ,
               16 , 17 , T(18)*c2 , T(19)*c2 , 20 , 21 , T(22)*c2 , T(23)*c2 ,
               24 , 25 , T(26)*c1 , T(27)*c1 , 28 , 29 , T(30)*c1 , T(31)*c1 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { crotz31.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    crotz62.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    std::vector< int > i8 = {  34 ,  35 ,  38 ,  39 ,  42 ,  43 ,  46 ,  47 ,
                               50 ,  51 ,  54 ,  55 ,  58 ,  59 ,  62 ,  63 ,
                               98 ,  99 , 102 , 103 , 106 , 107 , 110 , 111 ,
                              114 , 115 , 118 , 119 , 122 , 123 , 126 , 127 ,
                              162 , 163 , 166 , 167 , 170 , 171 , 174 , 175 ,
                              178 , 179 , 182 , 183 , 186 , 187 , 190 , 191 ,
                              226 , 227 , 230 , 231 , 234 , 235 , 238 , 239 ,
                              242 , 243 , 246 , 247 , 250 , 251 , 254 , 255 } ;
    V check8 = v8 ; for ( auto& i: i8 ) { check8[i-32] *= c1 ; check8[i] *= c2;}
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { crotz62.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    crotz62.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    check8 = v8 ; for ( auto& i: i8 ) { check8[i-32] *= c2 ; check8[i] *= c1 ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { crotz62.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

  //
  // CNOTs |0> controlled
  //

  // control < target
  {
    qclab::qgates::CRotationZ< T >  crotz01( 0 , 1 , pi/2 , 0 ) ;
    qclab::qgates::CRotationZ< T >  crotz12( 1 , 2 , pi/2 , 0 ) ;
    qclab::qgates::CRotationZ< T >  crotz02( 0 , 2 , pi/2 , 0 ) ;
    qclab::qgates::CRotationZ< T >  crotz13( 1 , 3 , pi/2 , 0 ) ;
    qclab::qgates::CRotationZ< T >  crotz26( 2 , 6 , pi/2 , 0 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    crotz01.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { T(3)*c1 , T(5)*c2 , 2 , 7 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { crotz01.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    crotz01.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check2 = { T(3)*c2 , T(5)*c1 , 2 , 7 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { crotz01.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    crotz01.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { T(3)*c1 , T(5)*c1 , T(2)*c2 , T(7)*c2 , 4 , 1 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz01.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz01.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { T(3)*c2 , T(5)*c2 , T(2)*c1 , T(7)*c1 , 4 , 1 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz01.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz12.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T(3)*c1 , T(5)*c2 , 2 , 7 , T(4)*c1 , T(1)*c2 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz12.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz12.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { T(3)*c2 , T(5)*c1 , 2 , 7 , T(4)*c2 , T(1)*c1 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz12.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz02.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T(3)*c1 , T(5)*c2 , T(2)*c1 , T(7)*c2 , 4 , 1 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz02.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz02.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { T(3)*c2 , T(5)*c1 , T(2)*c2 , T(7)*c1 , 4 , 1 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz02.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    crotz12.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { T(3)*c1 , T(5)*c1 , T(2)*c2 , T(7)*c2 , 4 , 1 , 8 , 3 ,
                 T(7)*c1 , T(2)*c1 , T(5)*c2 , T(6)*c2 , 8 , 9 , 5 , 1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { crotz12.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    crotz12.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { T(3)*c2 , T(5)*c2 , T(2)*c1 , T(7)*c1 , 4 , 1 , 8 , 3 ,
               T(7)*c2 , T(2)*c2 , T(5)*c1 , T(6)*c1 , 8 , 9 , 5 , 1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { crotz12.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    crotz13.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = { T( 0)*c1 , T( 1)*c1 , T( 2)*c2 , T( 3)*c2 ,
                 T( 4)*c1 , T( 5)*c1 , T( 6)*c2 , T( 7)*c2 ,
                  8 ,  9 , 10 , 11 , 12 , 13 , 14 , 15 ,
                 T(16)*c1 , T(17)*c1 , T(18)*c2 , T(19)*c2 ,
                 T(20)*c1 , T(21)*c1 , T(22)*c2 , T(23)*c2 ,
                 24 , 25 , 26 , 27 , 28 , 29 , 30 , 31 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { crotz13.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    crotz13.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    check5 = { T( 0)*c2 , T( 1)*c2 , T( 2)*c1 , T( 3)*c1 ,
               T( 4)*c2 , T( 5)*c2 , T( 6)*c1 , T( 7)*c1 ,
                8 ,  9 , 10 , 11 , 12 , 13 , 14 , 15 ,
               T(16)*c2 , T(17)*c2 , T(18)*c1 , T(19)*c1 ,
               T(20)*c2 , T(21)*c2 , T(22)*c1 , T(23)*c1 ,
               24 , 25 , 26 , 27 , 28 , 29 , 30 , 31 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { crotz13.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    crotz26.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    std::vector< int > i8 = {   2 ,   3 ,   6 ,   7 ,  10 ,  11 ,  14 ,  15 ,
                               18 ,  19 ,  22 ,  23 ,  26 ,  27 ,  30 ,  31 ,
                               66 ,  67 ,  70 ,  71 ,  74 ,  75 ,  78 ,  79 ,
                               82 ,  83 ,  86 ,  87 ,  90 ,  91 ,  94 ,  95 ,
                              130 , 131 , 134 , 135 , 138 , 139 , 142 , 143 ,
                              146 , 147 , 150 , 151 , 154 , 155 , 158 , 159 ,
                              194 , 195 , 198 , 199 , 202 , 203 , 206 , 207 ,
                              210 , 211 , 214 , 215 , 218 , 219 , 222 , 223 } ;
    V check8 = v8 ; for ( auto& i: i8 ) { check8[i-2] *= c1 ; check8[i] *= c2 ;}
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { crotz26.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    crotz26.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    check8 = v8 ; for ( auto& i: i8 ) { check8[i-2] *= c2 ; check8[i] *= c1 ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { crotz26.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

  // control > target
  {
    qclab::qgates::CRotationZ< T >  crotz10( 1 , 0 , pi/2 , 0 ) ;
    qclab::qgates::CRotationZ< T >  crotz21( 2 , 1 , pi/2 , 0 ) ;
    qclab::qgates::CRotationZ< T >  crotz20( 2 , 0 , pi/2 , 0 ) ;
    qclab::qgates::CRotationZ< T >  crotz31( 3 , 1 , pi/2 , 0 ) ;
    qclab::qgates::CRotationZ< T >  crotz62( 6 , 2 , pi/2 , 0 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    crotz10.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { T(3)*c1 , 5 , T(2)*c2 , 7 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { crotz10.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    crotz10.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check2 = { T(3)*c2 , 5 , T(2)*c1 , 7 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { crotz10.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    crotz10.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { T(3)*c1 , T(5)*c1 , 2 , 7 , T(4)*c2 , T(1)*c2 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz10.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz10.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { T(3)*c2 , T(5)*c2 , 2 , 7 , T(4)*c1 , T(1)*c1 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz10.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz21.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T(3)*c1 , 5 , T(2)*c2 , 7 , T(4)*c1 , 1 , T(8)*c2 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz21.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz21.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { T(3)*c2 , 5 , T(2)*c1 , 7 , T(4)*c2 , 1 , T(8)*c1 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz21.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz20.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T(3)*c1 , 5 , T(2)*c1 , 7 , T(4)*c2 , 1 , T(8)*c2 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz20.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    crotz20.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { T(3)*c2 , 5 , T(2)*c2 , 7 , T(4)*c1 , 1 , T(8)*c1 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { crotz20.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    crotz21.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { T(3)*c1 , T(5)*c1 , 2 , 7 , T(4)*c2 , T(1)*c2 , 8 , 3 ,
                 T(7)*c1 , T(2)*c1 , 5 , 6 , T(8)*c2 , T(9)*c2 , 5 , 1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { crotz21.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    crotz21.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { T(3)*c2 , T(5)*c2 , 2 , 7 , T(4)*c1 , T(1)*c1 , 8 , 3 ,
               T(7)*c2 , T(2)*c2 , 5 , 6 , T(8)*c1 , T(9)*c1 , 5 , 1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { crotz21.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    crotz31.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = { T( 0)*c1 , T( 1)*c1 ,  2 ,  3 , T( 4)*c1 , T( 5)*c1 ,  6 ,  7 ,
                 T( 8)*c2 , T( 9)*c2 , 10 , 11 , T(12)*c2 , T(13)*c2 , 14 , 15 ,
                 T(16)*c1 , T(17)*c1 , 18 , 19 , T(20)*c1 , T(21)*c1 , 22 , 23 ,
                 T(24)*c2 , T(25)*c2 , 26 , 27 , T(28)*c2 , T(29)*c2 , 30 , 31};
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { crotz31.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    crotz31.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    check5 = { T( 0)*c2 , T( 1)*c2 ,  2 ,  3 , T( 4)*c2 , T( 5)*c2 ,  6 ,  7 ,
               T( 8)*c1 , T( 9)*c1 , 10 , 11 , T(12)*c1 , T(13)*c1 , 14 , 15 ,
               T(16)*c2 , T(17)*c2 , 18 , 19 , T(20)*c2 , T(21)*c2 , 22 , 23 ,
               T(24)*c1 , T(25)*c1 , 26 , 27 , T(28)*c1 , T(29)*c1 , 30 , 31 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { crotz31.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    crotz62.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    std::vector< int > i8 = {  32 ,  33 ,  36 ,  37 ,  40 ,  41 ,  44 ,  45 ,
                               48 ,  49 ,  52 ,  53 ,  56 ,  57 ,  60 ,  61 ,
                               96 ,  97 , 100 , 101 , 104 , 105 , 108 , 109 ,
                              112 , 113 , 116 , 117 , 120 , 121 , 124 , 125 ,
                              160 , 161 , 164 , 165 , 168 , 169 , 172 , 173 ,
                              176 , 177 , 180 , 181 , 184 , 185 , 188 , 189 ,
                              224 , 225 , 228 , 229 , 232 , 233 , 236 , 237 ,
                              240 , 241 , 244 , 245 , 248 , 249 , 252 , 253 } ;
    V check8 = v8 ; for ( auto& i: i8 ) { check8[i-32] *= c1 ; check8[i] *= c2;}
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { crotz62.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    crotz62.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    check8 = v8 ; for ( auto& i: i8 ) { check8[i-32] *= c2 ; check8[i] *= c1 ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { crotz62.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
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

