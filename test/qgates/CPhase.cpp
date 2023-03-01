#include <gtest/gtest.h>
#include "qclab/qgates/CPhase.hpp"
#include <numeric>

template <typename T>
void test_qclab_qgates_CPhase() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R tol = 10 * std::numeric_limits< R >::epsilon() ;

  //
  // CPhases |1> controlled
  //
  {
    qclab::qgates::CPhase< T >  cphase ;

    EXPECT_EQ( cphase.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( cphase.fixed() ) ;          // fixed
    EXPECT_TRUE( cphase.controlled() ) ;      // controlled
    EXPECT_EQ( cphase.control() , 0 ) ;       // control
    EXPECT_EQ( cphase.target() , 1 ) ;        // target
    EXPECT_EQ( cphase.controlState() , 1 ) ;  // controlState

    // matrix
    auto eye = qclab::dense::eye< T >( 4 ) ;
    EXPECT_TRUE( cphase.matrix() == eye ) ;

    // qubit
    EXPECT_EQ( cphase.qubit() , 0 ) ;

    // qubits
    auto qubits = cphase.qubits() ;
    EXPECT_EQ( qubits.size() , 2 ) ;
    EXPECT_EQ( qubits[0] , 0 ) ;
    EXPECT_EQ( qubits[1] , 1 ) ;
    int qnew[] = { 5 , 3 } ;
    cphase.setQubits( &qnew[0] ) ;
    EXPECT_EQ( cphase.qubits()[0] , 3 ) ;
    EXPECT_EQ( cphase.qubits()[1] , 5 ) ;
    qnew[0] = 0 ;
    qnew[1] = 1 ;
    cphase.setQubits( &qnew[0] ) ;

    // print
    cphase.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cphase.toQASM( qasm ) , 0 ) ;
    auto str_theta = qclab::qasm( cphase.theta() ) ;
    std::string qasm_check = "cp(" + str_theta + ") q[0], q[1];\n" ;
    EXPECT_EQ( qasm.str() , qasm_check ) ;
    std::cout << qasm.str() ;

    // gate
    qclab::qgates::Phase< T >  P ;
    EXPECT_TRUE( *cphase.gate() == P ) ;
    EXPECT_TRUE( cphase.gate()->matrix() == P.matrix() ) ;

    // operators == and !=
    EXPECT_TRUE(  cphase != P ) ;
    EXPECT_FALSE( cphase == P ) ;
    qclab::qgates::CPhase< T >  cphase2 ;
    EXPECT_TRUE(  cphase == cphase2 ) ;
    EXPECT_FALSE( cphase != cphase2 ) ;

    // setControl, setTarget, setControlState
    cphase.setControl( 3 ) ;
    EXPECT_EQ( cphase.control() , 3 ) ;
    cphase.setTarget( 5 ) ;
    EXPECT_EQ( cphase.target() , 5 ) ;
    EXPECT_TRUE(  cphase == cphase2 ) ;
    EXPECT_FALSE( cphase != cphase2 ) ;

    cphase.setControl( 4 ) ;
    EXPECT_EQ( cphase.control() , 4 ) ;
    cphase.setTarget( 1 ) ;
    EXPECT_EQ( cphase.target() , 1 ) ;
    EXPECT_TRUE(  cphase != cphase2 ) ;
    EXPECT_FALSE( cphase == cphase2 ) ;

    qubits[0] = 1 ;
    qubits[1] = 2 ;
    cphase.setQubits( &qubits[0] ) ;
    EXPECT_EQ( cphase.control() , 1 ) ;
    EXPECT_EQ( cphase.target() , 2 ) ;
    EXPECT_TRUE(  cphase == cphase2 ) ;
    EXPECT_FALSE( cphase != cphase2 ) ;

    cphase.setControl( 0 ) ;
    EXPECT_EQ( cphase.control() , 0 ) ;
    cphase.setTarget( 1 ) ;
    EXPECT_EQ( cphase.target() , 1 ) ;
    cphase.setControlState( 0 ) ;
    EXPECT_EQ( cphase.controlState() , 0 ) ;
    EXPECT_TRUE(  cphase != cphase2 ) ;
    EXPECT_FALSE( cphase == cphase2 ) ;

    // makeFixed, makeVariable
    cphase.makeFixed() ;
    EXPECT_TRUE( cphase.fixed() ) ;
    cphase.makeVariable() ;
    EXPECT_FALSE( cphase.fixed() ) ;

    // angle, theta, sin, cos
    qclab::QAngle< R > angle ;
    EXPECT_TRUE( cphase.angle() == angle ) ;
    EXPECT_EQ( cphase.theta() , 0 ) ;
    EXPECT_EQ( cphase.cos() , 1 ) ;
    EXPECT_EQ( cphase.sin() , 0 ) ;

    // update(angle)
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    angle.update( pi/4 ) ;
    cphase.update( angle ) ;
    EXPECT_NEAR( cphase.theta() , pi/4 , tol ) ;
    EXPECT_NEAR( cphase.cos() , cos , tol ) ;
    EXPECT_NEAR( cphase.sin() , sin , tol ) ;

    // update(theta)
    cphase.update( pi/2 ) ;
    EXPECT_NEAR( cphase.theta() , pi/2 , tol ) ;
    EXPECT_NEAR( cphase.cos() , 0 , tol ) ;
    EXPECT_NEAR( cphase.sin() , 1 , tol ) ;

    // update(cos,sin)
    cphase.update( cos , sin ) ;
    EXPECT_NEAR( cphase.theta() , pi/4 , tol ) ;
    EXPECT_NEAR( cphase.cos() , cos , tol ) ;
    EXPECT_NEAR( cphase.sin() , sin , tol ) ;
  }

  {
    const R theta = pi/4 ;
    qclab::QAngle< R >  angle( theta ) ;
    qclab::qgates::CPhase< T >  cphase( 1 , 3 , angle ) ;

    EXPECT_EQ( cphase.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( cphase.fixed() ) ;          // fixed
    EXPECT_TRUE( cphase.controlled() ) ;      // controlled
    EXPECT_EQ( cphase.control() , 1 ) ;       // control
    EXPECT_EQ( cphase.target() , 3 ) ;        // target
    EXPECT_EQ( cphase.controlState() , 1 ) ;  // controlState

    EXPECT_NEAR( cphase.theta() , theta , tol ) ;
    EXPECT_NEAR( cphase.cos() , std::cos( theta ) , tol ) ;
    EXPECT_NEAR( cphase.sin() , std::sin( theta ) , tol ) ;
  }

  {
    const R theta = pi/4 ;
    qclab::qgates::CPhase< T >  cphase( 1 , 3 , theta ) ;

    EXPECT_EQ( cphase.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( cphase.fixed() ) ;          // fixed
    EXPECT_TRUE( cphase.controlled() ) ;      // controlled
    EXPECT_EQ( cphase.control() , 1 ) ;       // control
    EXPECT_EQ( cphase.target() , 3 ) ;        // target
    EXPECT_EQ( cphase.controlState() , 1 ) ;  // controlState

    EXPECT_NEAR( cphase.theta() , theta , tol ) ;
    EXPECT_NEAR( cphase.cos() , std::cos( theta ) , tol ) ;
    EXPECT_NEAR( cphase.sin() , std::sin( theta ) , tol ) ;
  }

  {
    const R theta = pi/4 ;
    const R cos = std::cos( theta ) ;
    const R sin = std::sin( theta ) ;
    qclab::qgates::CPhase< T >  cphase( 1 , 3 , cos , sin ) ;

    EXPECT_EQ( cphase.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_FALSE( cphase.fixed() ) ;          // fixed
    EXPECT_TRUE( cphase.controlled() ) ;      // controlled
    EXPECT_EQ( cphase.control() , 1 ) ;       // control
    EXPECT_EQ( cphase.target() , 3 ) ;        // target
    EXPECT_EQ( cphase.controlState() , 1 ) ;  // controlState

    EXPECT_NEAR( cphase.theta() , theta , tol ) ;
    EXPECT_NEAR( cphase.cos() , cos , tol ) ;
    EXPECT_NEAR( cphase.sin() , sin , tol ) ;
  }

  using V = std::vector< T > ;

  const T c  = T( std::cos( pi/4 ) ,  std::sin( pi/4 ) ) ;
  const T cc = T( std::cos( pi/4 ) , -std::sin( pi/4 ) ) ;

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
    qclab::qgates::CPhase< T >  cphase01( 0 , 1 , pi/4 ) ;
    qclab::qgates::CPhase< T >  cphase12( 1 , 2 , pi/4 ) ;
    qclab::qgates::CPhase< T >  cphase02( 0 , 2 , pi/4 ) ;
    qclab::qgates::CPhase< T >  cphase13( 1 , 3 , pi/4 ) ;
    qclab::qgates::CPhase< T >  cphase26( 2 , 6 , pi/4 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cphase01.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 3 , 5 , 2 , T(7)*c } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cphase01.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cphase01.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check2 = { 3 , 5 , 2 , T(7)*std::conj(c) } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cphase01.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    cphase01.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 3 , 5 , 2 , 7 , 4 , 1 , T(8)*c , T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase01.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase01.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , 7 , 4 , 1 , T(8)*std::conj(c) , T(3)*std::conj(c) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase01.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase12.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , T(7)*c , 4 , 1 , 8 , T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase12.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase12.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , T(7)*std::conj(c) , 4 , 1 , 8 , T(3)*std::conj(c) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase12.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase02.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , 7 , 4 , T(1)*c , 8 , T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase02.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase02.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , 7 , 4 , T(1)*std::conj(c) , 8 , T(3)*std::conj(c) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase02.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    cphase12.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { 3 , 5 , 2 , 7 , 4 , 1 , T(8)*c , T(3)*c ,
                 7 , 2 , 5 , 6 , 8 , 9 , T(5)*c , T(1)*c } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cphase12.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cphase12.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { 3 , 5 , 2 , 7 , 4 , 1 , T(8)*std::conj(c) , T(3)*std::conj(c) ,
               7 , 2 , 5 , 6 , 8 , 9 , T(5)*std::conj(c) , T(1)*std::conj(c) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cphase12.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    cphase13.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = {  0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,
                  8 ,  9 , T(10)*c , T(11)*c , 12 , 13 , T(14)*c , T(15)*c ,
                 16 , 17 , 18 , 19 , 20 , 21 , 22 , 23 ,
                 24 , 25 , T(26)*c , T(27)*c , 28 , 29 , T(30)*c , T(31)*c } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cphase13.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cphase13.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    check5 = {  0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,
                8 ,  9 , T(10)*cc , T(11)*cc , 12 , 13 , T(14)*cc , T(15)*cc ,
               16 , 17 , 18 , 19 , 20 , 21 , 22 , 23 ,
               24 , 25 , T(26)*cc , T(27)*cc , 28 , 29 , T(30)*cc , T(31)*cc } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cphase13.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    cphase26.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    std::vector< int > i8 = {  34 ,  35 ,  38 ,  39 ,  42 ,  43 ,  46 ,  47 ,
                               50 ,  51 ,  54 ,  55 ,  58 ,  59 ,  62 ,  63 ,
                               98 ,  99 , 102 , 103 , 106 , 107 , 110 , 111 ,
                              114 , 115 , 118 , 119 , 122 , 123 , 126 , 127 ,
                              162 , 163 , 166 , 167 , 170 , 171 , 174 , 175 ,
                              178 , 179 , 182 , 183 , 186 , 187 , 190 , 191 ,
                              226 , 227 , 230 , 231 , 234 , 235 , 238 , 239 ,
                              242 , 243 , 246 , 247 , 250 , 251 , 254 , 255 } ;
    V check8 = v8 ; for ( auto& i: i8 ) check8[i] *= c ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cphase26.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cphase26.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    check8 = v8 ; for ( auto& i: i8 ) check8[i] *= cc ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cphase26.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

  // control > target
  {
    qclab::qgates::CPhase< T >  cphase10( 1 , 0 , pi/4 ) ;
    qclab::qgates::CPhase< T >  cphase21( 2 , 1 , pi/4 ) ;
    qclab::qgates::CPhase< T >  cphase20( 2 , 0 , pi/4 ) ;
    qclab::qgates::CPhase< T >  cphase31( 3 , 1 , pi/4 ) ;
    qclab::qgates::CPhase< T >  cphase62( 6 , 2 , pi/4 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cphase10.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 3 , 5 , 2 , T(7)*c } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cphase10.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cphase10.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check2 = { 3 , 5 , 2 , T(7)*std::conj(c) } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cphase10.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    cphase10.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 3 , 5 , 2 , 7 , 4 , 1 , T(8)*c , T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase10.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase10.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , 7 , 4 , 1 , T(8)*std::conj(c) , T(3)*std::conj(c) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase10.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase21.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , T(7)*c , 4 , 1 , 8 , T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase21.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase21.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , T(7)*std::conj(c) , 4 , 1 , 8 , T(3)*std::conj(c) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase21.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase20.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , 7 , 4 , T(1)*c , 8 , T(3)*c } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase20.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase20.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , 7 , 4 , T(1)*std::conj(c) , 8 , T(3)*std::conj(c) } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase20.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    cphase21.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { 3 , 5 , 2 , 7 , 4 , 1 , T(8)*c , T(3)*c ,
                 7 , 2 , 5 , 6 , 8 , 9 , T(5)*c , T(1)*c } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cphase21.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cphase21.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { 3 , 5 , 2 , 7 , 4 , 1 , T(8)*std::conj(c) , T(3)*std::conj(c) ,
               7 , 2 , 5 , 6 , 8 , 9 , T(5)*std::conj(c) , T(1)*std::conj(c) } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cphase21.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    cphase31.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = {  0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,
                  8 ,  9 , T(10)*c , T(11)*c , 12 , 13 , T(14)*c , T(15)*c ,
                 16 , 17 , 18 , 19 , 20 , 21 , 22 , 23 ,
                 24 , 25 , T(26)*c , T(27)*c , 28 , 29 , T(30)*c , T(31)*c } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cphase31.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cphase31.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    check5 = {  0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,
                8 ,  9 , T(10)*cc , T(11)*cc , 12 , 13 , T(14)*cc , T(15)*cc ,
               16 , 17 , 18 , 19 , 20 , 21 , 22 , 23 ,
               24 , 25 , T(26)*cc , T(27)*cc , 28 , 29 , T(30)*cc , T(31)*cc } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cphase31.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    cphase62.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    std::vector< int > i8 = {  34 ,  35 ,  38 ,  39 ,  42 ,  43 ,  46 ,  47 ,
                               50 ,  51 ,  54 ,  55 ,  58 ,  59 ,  62 ,  63 ,
                               98 ,  99 , 102 , 103 , 106 , 107 , 110 , 111 ,
                              114 , 115 , 118 , 119 , 122 , 123 , 126 , 127 ,
                              162 , 163 , 166 , 167 , 170 , 171 , 174 , 175 ,
                              178 , 179 , 182 , 183 , 186 , 187 , 190 , 191 ,
                              226 , 227 , 230 , 231 , 234 , 235 , 238 , 239 ,
                              242 , 243 , 246 , 247 , 250 , 251 , 254 , 255 } ;
    V check8 = v8 ; for ( auto& i: i8 ) check8[i] *= c ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cphase62.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cphase62.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    check8 = v8 ; for ( auto& i: i8 ) check8[i] *= cc ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cphase62.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

  //
  // CNOTs |0> controlled
  //

  // control < target
  {
    qclab::qgates::CPhase< T >  cphase01( 0 , 1 , pi/4 , 0 ) ;
    qclab::qgates::CPhase< T >  cphase12( 1 , 2 , pi/4 , 0 ) ;
    qclab::qgates::CPhase< T >  cphase02( 0 , 2 , pi/4 , 0 ) ;
    qclab::qgates::CPhase< T >  cphase13( 1 , 3 , pi/4 , 0 ) ;
    qclab::qgates::CPhase< T >  cphase26( 2 , 6 , pi/4 , 0 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cphase01.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 3 , T(5)*c , 2 , 7 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cphase01.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cphase01.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check2 = { 3 , T(5)*std::conj(c) , 2 , 7 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cphase01.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    cphase01.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 3 , 5 , T(2)*c , T(7)*c , 4 , 1 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase01.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase01.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , T(2)*std::conj(c) , T(7)*std::conj(c) , 4 , 1 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase01.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase12.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , T(5)*c , 2 , 7 , 4 , T(1)*c , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase12.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase12.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , T(5)*std::conj(c) , 2 , 7 , 4 , T(1)*std::conj(c) , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase12.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase02.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , T(5)*c , 2 , T(7)*c , 4 , 1 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase02.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase02.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , T(5)*std::conj(c) , 2 , T(7)*std::conj(c) , 4 , 1 , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase02.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    cphase12.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { 3 , 5 , T(2)*c , T(7)*c , 4 , 1 , 8 , 3 ,
                 7 , 2 , T(5)*c , T(6)*c , 8 , 9 , 5 , 1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cphase12.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cphase12.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { 3 , 5 , T(2)*std::conj(c) , T(7)*std::conj(c) , 4 , 1 , 8 , 3 ,
               7 , 2 , T(5)*std::conj(c) , T(6)*std::conj(c) , 8 , 9 , 5 , 1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cphase12.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    cphase13.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = {  0 ,  1 , T( 2)*c , T( 3)*c ,  4 ,  5 , T( 6)*c , T( 7)*c ,
                  8 ,  9 , 10 , 11 , 12 , 13 , 14 , 15 ,
                 16 , 17 , T(18)*c , T(19)*c , 20 , 21 , T(22)*c , T(23)*c ,
                 24 , 25 , 26 , 27 , 28 , 29 , 30 , 31 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cphase13.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cphase13.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    check5 = {  0 ,  1 , T( 2)*cc , T( 3)*cc ,  4 ,  5 , T( 6)*cc , T( 7)*cc ,
                8 ,  9 , 10 , 11 , 12 , 13 , 14 , 15 ,
               16 , 17 , T(18)*cc , T(19)*cc , 20 , 21 , T(22)*cc , T(23)*cc ,
               24 , 25 , 26 , 27 , 28 , 29 , 30 , 31 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cphase13.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    cphase26.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    std::vector< int > i8 = {   2 ,   3 ,   6 ,   7 ,  10 ,  11 ,  14 ,  15 ,
                               18 ,  19 ,  22 ,  23 ,  26 ,  27 ,  30 ,  31 ,
                               66 ,  67 ,  70 ,  71 ,  74 ,  75 ,  78 ,  79 ,
                               82 ,  83 ,  86 ,  87 ,  90 ,  91 ,  94 ,  95 ,
                              130 , 131 , 134 , 135 , 138 , 139 , 142 , 143 ,
                              146 , 147 , 150 , 151 , 154 , 155 , 158 , 159 ,
                              194 , 195 , 198 , 199 , 202 , 203 , 206 , 207 ,
                              210 , 211 , 214 , 215 , 218 , 219 , 222 , 223 } ;
    V check8 = v8 ; for ( auto& i: i8 ) check8[i] *= c ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cphase26.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cphase26.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    check8 = v8 ; for ( auto& i: i8 ) check8[i] *= cc ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cphase26.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

  // control > target
  {
    qclab::qgates::CPhase< T >  cphase10( 1 , 0 , pi/4 , 0 ) ;
    qclab::qgates::CPhase< T >  cphase21( 2 , 1 , pi/4 , 0 ) ;
    qclab::qgates::CPhase< T >  cphase20( 2 , 0 , pi/4 , 0 ) ;
    qclab::qgates::CPhase< T >  cphase31( 3 , 1 , pi/4 , 0 ) ;
    qclab::qgates::CPhase< T >  cphase62( 6 , 2 , pi/4 , 0 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cphase10.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 3 , 5 , T(2)*c , 7 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cphase10.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cphase10.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check2 = { 3 , 5 , T(2)*std::conj(c) , 7 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cphase10.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    cphase10.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 3 , 5 , 2 , 7 , T(4)*c , T(1)*c , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase10.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase10.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , 7 , T(4)*std::conj(c) , T(1)*std::conj(c) , 8 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase10.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase21.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , T(2)*c , 7 , 4 , 1 , T(8)*c , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase21.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase21.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , T(2)*std::conj(c) , 7 , 4 , 1 , T(8)*std::conj(c) , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase21.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase20.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , 7 , T(4)*c , 1 , T(8)*c , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase20.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cphase20.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { 3 , 5 , 2 , 7 , T(4)*std::conj(c) , 1 , T(8)*std::conj(c) , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cphase20.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    cphase21.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { 3 , 5 , 2 , 7 , T(4)*c , T(1)*c , 8 , 3 ,
                 7 , 2 , 5 , 6 , T(8)*c , T(9)*c , 5 , 1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cphase21.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cphase21.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { 3 , 5 , 2 , 7 , T(4)*std::conj(c) , T(1)*std::conj(c) , 8 , 3 ,
               7 , 2 , 5 , 6 , T(8)*std::conj(c) , T(9)*std::conj(c) , 5 , 1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cphase21.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    cphase31.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = {  0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,
                 T( 8)*c , T( 9)*c , 10 , 11 , T(12)*c , T(13)*c , 14 , 15 ,
                 16 , 17 , 18 , 19 , 20 , 21 , 22 , 23 ,
                 T(24)*c , T(25)*c , 26 , 27 , T(28)*c , T(29)*c , 30 , 31 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cphase31.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cphase31.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    check5 = {  0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,
               T( 8)*cc , T( 9)*cc , 10 , 11 , T(12)*cc , T(13)*cc , 14 , 15 ,
               16 , 17 , 18 , 19 , 20 , 21 , 22 , 23 ,
               T(24)*cc , T(25)*cc , 26 , 27 , T(28)*cc , T(29)*cc , 30 , 31 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cphase31.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    cphase62.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    std::vector< int > i8 = {  32 ,  33 ,  36 ,  37 ,  40 ,  41 ,  44 ,  45 ,
                               48 ,  49 ,  52 ,  53 ,  56 ,  57 ,  60 ,  61 ,
                               96 ,  97 , 100 , 101 , 104 , 105 , 108 , 109 ,
                              112 , 113 , 116 , 117 , 120 , 121 , 124 , 125 ,
                              160 , 161 , 164 , 165 , 168 , 169 , 172 , 173 ,
                              176 , 177 , 180 , 181 , 184 , 185 , 188 , 189 ,
                              224 , 225 , 228 , 229 , 232 , 233 , 236 , 237 ,
                              240 , 241 , 244 , 245 , 248 , 249 , 252 , 253 } ;
    V check8 = v8 ; for ( auto& i: i8 ) check8[i] *= c ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cphase62.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cphase62.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    check8 = v8 ; for ( auto& i: i8 ) check8[i] *= cc ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cphase62.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

}


/*
 * complex float
 */
TEST( qclab_qgates_CPhase , complex_float ) {
  test_qclab_qgates_CPhase< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_CPhase , complex_double ) {
  test_qclab_qgates_CPhase< std::complex< double > >() ;
}

