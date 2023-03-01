#include <gtest/gtest.h>
#include "qclab/qgates/RotationY.hpp"

template <typename T>
bool check_qclab_qgates_RotationY( const std::vector< T >& v1 ,
                                   const std::vector< T >& v2 ) {

  using R = qclab::real_t< T > ;
  const R tol = 100 * std::numeric_limits< R >::epsilon() ;

  assert( v1.size() == v2.size() ) ;
  for ( size_t i = 0; i < v1.size(); ++i ) {
    if ( std::abs( v1[i] - v2[i] ) > tol ) return false ;
  }
  return true ;

}

template <typename T>
void test_qclab_qgates_RotationY() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R tol = 10 * std::numeric_limits< R >::epsilon() ;

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

    // update(rot)
    qclab::QRotation< R >  new_rot( 1 ) ;
    Ry.update( new_rot ) ;
    EXPECT_NEAR( Ry.theta() , 1 , tol ) ;
    EXPECT_NEAR( Ry.cos() , std::cos( 0.5 ) , tol ) ;
    EXPECT_NEAR( Ry.sin() , std::sin( 0.5 ) , tol ) ;

    // update(theta)
    Ry.update( pi/2 ) ;
    EXPECT_NEAR( Ry.theta() , pi/2 , tol ) ;
    EXPECT_NEAR( Ry.cos() , std::cos( pi/4 ) , tol ) ;
    EXPECT_NEAR( Ry.sin() , std::sin( pi/4 ) , tol ) ;

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
    EXPECT_NEAR( Ry.theta() , 2*(pi/3) , tol ) ;

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

    EXPECT_NEAR( Ry.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Ry.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Ry.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    qclab::qgates::RotationY< T >  Ry( 5 , pi/2 ) ;

    EXPECT_EQ( Ry.nbQubits() , 1 ) ;                    // nbQubits
    EXPECT_FALSE( Ry.fixed() ) ;                        // fixed
    EXPECT_FALSE( Ry.controlled() ) ;                   // controlled
    EXPECT_EQ( Ry.qubit() , 5 ) ;                       // qubit

    EXPECT_NEAR( Ry.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Ry.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Ry.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    qclab::qgates::RotationY< T >  Ry( 5 , pi/2 , true ) ;

    EXPECT_EQ( Ry.nbQubits() , 1 ) ;                    // nbQubits
    EXPECT_TRUE( Ry.fixed() ) ;                         // fixed
    EXPECT_FALSE( Ry.controlled() ) ;                   // controlled
    EXPECT_EQ( Ry.qubit() , 5 ) ;                       // qubit

    EXPECT_NEAR( Ry.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Ry.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Ry.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationY< T >  Ry( cos , sin ) ;

    EXPECT_EQ( Ry.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_FALSE( Ry.fixed() ) ;              // fixed
    EXPECT_FALSE( Ry.controlled() ) ;         // controlled
    EXPECT_EQ( Ry.qubit() , 0 ) ;             // qubit

    EXPECT_NEAR( Ry.theta() , pi/2 , tol ) ;  // theta
    EXPECT_NEAR( Ry.cos() , cos , tol ) ;     // cos
    EXPECT_NEAR( Ry.sin() , sin , tol ) ;     // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationY< T >  Ry( 5 , cos , sin ) ;

    EXPECT_EQ( Ry.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_FALSE( Ry.fixed() ) ;              // fixed
    EXPECT_FALSE( Ry.controlled() ) ;         // controlled
    EXPECT_EQ( Ry.qubit() , 5 ) ;             // qubit

    EXPECT_NEAR( Ry.theta() , pi/2 , tol ) ;  // theta
    EXPECT_NEAR( Ry.cos() , cos , tol ) ;     // cos
    EXPECT_NEAR( Ry.sin() , sin , tol ) ;     // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationY< T >  Ry( 5 , cos , sin , true ) ;

    EXPECT_EQ( Ry.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_TRUE( Ry.fixed() ) ;               // fixed
    EXPECT_FALSE( Ry.controlled() ) ;         // controlled
    EXPECT_EQ( Ry.qubit() , 5 ) ;             // qubit

    EXPECT_NEAR( Ry.theta() , pi/2 , tol ) ;  // theta
    EXPECT_NEAR( Ry.cos() , cos , tol ) ;     // cos
    EXPECT_NEAR( Ry.sin() , sin , tol ) ;     // sin
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >       angle1( theta1 ) ;
    const qclab::QRotation< R >    rot1( angle1 ) ;
    qclab::qgates::RotationY< T >  R1( rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >       angle2( theta2 ) ;
    const qclab::QRotation< R >    rot2( angle2 ) ;
    qclab::qgates::RotationY< T >  R2( rot2 ) ;

    // operator *=
    R1 *= R2 ;
    qclab::QAngle< R >  angle = angle1 + angle2 ;
    R cos   = angle.cos() ;
    R sin   = angle.sin() ;
    R theta = angle.theta() ;
    EXPECT_NEAR( R1.cos() , cos , tol ) ;                 // cos
    EXPECT_NEAR( R1.sin() , sin , tol ) ;                 // sin
    EXPECT_NEAR( R1.theta() , 2*theta , tol ) ;           // theta
    EXPECT_NEAR( R2.cos() , std::cos( theta2 ) , tol ) ;  // cos
    EXPECT_NEAR( R2.sin() , std::sin( theta2 ) , tol ) ;  // sin
    EXPECT_NEAR( R2.theta() , 2*theta2 , tol ) ;          // theta
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >       angle1( theta1 ) ;
    const qclab::QRotation< R >    rot1( angle1 ) ;
    qclab::qgates::RotationY< T >  R1( rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >       angle2( theta2 ) ;
    const qclab::QRotation< R >    rot2( angle2 ) ;
    qclab::qgates::RotationY< T >  R2( rot2 ) ;

    // operator /=
    R1 /= R2 ;
    qclab::QAngle< R >  angle = angle1 - angle2 ;
    R cos   = angle.cos() ;
    R sin   = angle.sin() ;
    R theta = angle.theta() ;
    EXPECT_NEAR( R1.cos() , cos , tol ) ;                 // cos
    EXPECT_NEAR( R1.sin() , sin , tol ) ;                 // sin
    EXPECT_NEAR( R1.theta() , 2*theta , tol ) ;           // theta
    EXPECT_NEAR( R2.cos() , std::cos( theta2 ) , tol ) ;  // cos
    EXPECT_NEAR( R2.sin() , std::sin( theta2 ) , tol ) ;  // sin
    EXPECT_NEAR( R2.theta() , 2*theta2 , tol ) ;          // theta
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >       angle1( theta1 ) ;
    const qclab::QRotation< R >    rot1( angle1 ) ;
    qclab::qgates::RotationY< T >  R1( rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >       angle2( theta2 ) ;
    const qclab::QRotation< R >    rot2( angle2 ) ;
    qclab::qgates::RotationY< T >  R2( rot2 ) ;

    // operator *
    qclab::qgates::RotationY< T > R12 = R1 * R2 ;
    qclab::QAngle< R >  angle = angle1 + angle2 ;
    R cos   = angle.cos() ;
    R sin   = angle.sin() ;
    R theta = angle.theta() ;
    EXPECT_NEAR( R12.cos() , cos , tol ) ;        // cos
    EXPECT_NEAR( R12.sin() , sin , tol ) ;        // sin
    EXPECT_NEAR( R12.theta() , 2*theta , tol ) ;  // theta
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >       angle1( theta1 ) ;
    const qclab::QRotation< R >    rot1( angle1 ) ;
    qclab::qgates::RotationY< T >  R1( rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >       angle2( theta2 ) ;
    const qclab::QRotation< R >    rot2( angle2 ) ;
    qclab::qgates::RotationY< T >  R2( rot2 ) ;

    // operator /
    qclab::qgates::RotationY< T > R12 = R1 / R2 ;
    qclab::QAngle< R >  angle = angle1 - angle2 ;
    R cos   = angle.cos() ;
    R sin   = angle.sin() ;
    R theta = angle.theta() ;
    EXPECT_NEAR( R12.cos() , cos , tol ) ;        // cos
    EXPECT_NEAR( R12.sin() , sin , tol ) ;        // sin
    EXPECT_NEAR( R12.theta() , 2*theta , tol ) ;  // theta
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >       angle1( theta1 ) ;
    const qclab::QRotation< R >    rot1( angle1 ) ;
    qclab::qgates::RotationY< T >  R1( rot1 ) ;
    R cos   = R1.cos() ;
    R sin   = R1.sin() ;
    R theta = R1.theta() ;

    // inv
    qclab::qgates::RotationY< T >  R2 = R1.inv() ;
    EXPECT_NEAR( R1.cos() ,  cos , tol ) ;      // cos
    EXPECT_NEAR( R2.cos() ,  cos , tol ) ;      // cos
    EXPECT_NEAR( R1.sin() ,  sin , tol ) ;      // sin
    EXPECT_NEAR( R2.sin() , -sin , tol ) ;      // sin
    EXPECT_NEAR( R1.theta() ,  theta , tol ) ;  // theta
    EXPECT_NEAR( R2.theta() , -theta , tol ) ;  // theta
  }

  // apply
  {
    using V = std::vector< T > ;
    const V v1 = { 1 , 2 } ;
    const V v2 = { 1 , 2 , 3 , 4 } ;
    const V v3 = { 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 } ;
    const V v4 = {  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  8 ,
                    9 , 10 , 11 , 12 , 13 , 14 , 15 , 16 } ;

    const R c = std::cos( pi/8 ) ;
    const R s = std::sin( pi/8 ) ;

    qclab::qgates::RotationY< T >  RY0( 0 , pi/4 ) ;
    qclab::qgates::RotationY< T >  RY1( 1 , pi/4 ) ;
    qclab::qgates::RotationY< T >  RY2( 2 , pi/4 ) ;
    qclab::qgates::RotationY< T >  RY3( 3 , pi/4 ) ;

    // nbQubits = 1
    auto vec1 = v1 ;
    RY0.apply( qclab::Op::NoTrans , 1 , vec1 ) ;
    V check1 = { T(1)*c - T(2)*s , T(1)*s + T(2)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec1 , check1 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ; T* vec1_ = vec1.data() ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { RY0.apply_device( qclab::Op::NoTrans , 1 , vec1_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec1 , check1 ) ) ;
  #endif

    vec1 = v1 ;
    RY0.apply( qclab::Op::Trans , 1 , vec1 ) ;
    check1 = { T(1)*c + T(2)*s , -T(1)*s + T(2)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec1 , check1 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { RY0.apply_device( qclab::Op::Trans , 1 , vec1_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec1 , check1 ) ) ;
  #endif

    vec1 = v1 ;
    RY0.apply( qclab::Op::ConjTrans , 1 , vec1 ) ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec1 , check1 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { RY0.apply_device( qclab::Op::ConjTrans , 1 , vec1_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec1 , check1 ) ) ;
  #endif

    // nbQubits = 2
    auto vec2 = v2 ;
    RY0.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { T(1)*c - T(3)*s , T(2)*c - T(4)*s ,
                 T(1)*s + T(3)*c , T(2)*s + T(4)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec2 , check2 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { RY0.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec2 , check2 ) ) ;
  #endif

    vec2 = v2 ;
    RY0.apply( qclab::Op::Trans , 2 , vec2 ) ;
    check2 = {  T(1)*c + T(3)*s ,  T(2)*c + T(4)*s ,
               -T(1)*s + T(3)*c , -T(2)*s + T(4)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec2 , check2 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { RY0.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec2 , check2 ) ) ;
  #endif

    vec2 = v2 ;
    RY0.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec2 , check2 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { RY0.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec2 , check2 ) ) ;
  #endif

    vec2 = v2 ;
    RY1.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    check2 = { T(1)*c - T(2)*s , T(1)*s + T(2)*c ,
               T(3)*c - T(4)*s , T(3)*s + T(4)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec2 , check2 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { RY1.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec2 , check2 ) ) ;
  #endif

    vec2 = v2 ;
    RY1.apply( qclab::Op::Trans , 2 , vec2 ) ;
    check2 = { T(1)*c + T(2)*s , -T(1)*s + T(2)*c ,
               T(3)*c + T(4)*s , -T(3)*s + T(4)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec2 , check2 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { RY1.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec2 , check2 ) ) ;
  #endif

    vec2 = v2 ;
    RY1.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec2 , check2 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { RY1.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec2 , check2 ) ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    RY0.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { T(1)*c - T(5)*s , T(2)*c - T(6)*s ,
                 T(3)*c - T(7)*s , T(4)*c - T(8)*s ,
                 T(1)*s + T(5)*c , T(2)*s + T(6)*c ,
                 T(3)*s + T(7)*c , T(4)*s + T(8)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RY0.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #endif

    vec3 = v3 ;
    RY0.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = {  T(1)*c + T(5)*s ,  T(2)*c + T(6)*s ,
                T(3)*c + T(7)*s ,  T(4)*c + T(8)*s ,
               -T(1)*s + T(5)*c , -T(2)*s + T(6)*c ,
               -T(3)*s + T(7)*c , -T(4)*s + T(8)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RY0.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #endif

    vec3 = v3 ;
    RY0.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RY0.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #endif

    vec3 = v3 ;
    RY1.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T(1)*c - T(3)*s , T(2)*c - T(4)*s ,
               T(1)*s + T(3)*c , T(2)*s + T(4)*c ,
               T(5)*c - T(7)*s , T(6)*c - T(8)*s ,
               T(5)*s + T(7)*c , T(6)*s + T(8)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RY1.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #endif

    vec3 = v3 ;
    RY1.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = {  T(1)*c + T(3)*s ,  T(2)*c + T(4)*s ,
               -T(1)*s + T(3)*c , -T(2)*s + T(4)*c ,
                T(5)*c + T(7)*s ,  T(6)*c + T(8)*s ,
               -T(5)*s + T(7)*c , -T(6)*s + T(8)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RY1.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #endif

    vec3 = v3 ;
    RY1.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RY1.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #endif

    vec3 = v3 ;
    RY2.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T(1)*c - T(2)*s , T(1)*s + T(2)*c ,
               T(3)*c - T(4)*s , T(3)*s + T(4)*c ,
               T(5)*c - T(6)*s , T(5)*s + T(6)*c ,
               T(7)*c - T(8)*s , T(7)*s + T(8)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RY2.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #endif

    vec3 = v3 ;
    RY2.apply( qclab::Op::Trans , 3 , vec3 ) ;
    check3 = { T(1)*c + T(2)*s , -T(1)*s + T(2)*c ,
               T(3)*c + T(4)*s , -T(3)*s + T(4)*c ,
               T(5)*c + T(6)*s , -T(5)*s + T(6)*c ,
               T(7)*c + T(8)*s , -T(7)*s + T(8)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RY2.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #endif

    vec3 = v3 ;
    RY2.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RY2.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec3 , check3 ) ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    RY0.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { T( 1)*c - T( 9)*s , T( 2)*c - T(10)*s ,
                 T( 3)*c - T(11)*s , T( 4)*c - T(12)*s ,
                 T( 5)*c - T(13)*s , T( 6)*c - T(14)*s ,
                 T( 7)*c - T(15)*s , T( 8)*c - T(16)*s ,
                 T( 1)*s + T( 9)*c , T( 2)*s + T(10)*c ,
                 T( 3)*s + T(11)*c , T( 4)*s + T(12)*c ,
                 T( 5)*s + T(13)*c , T( 6)*s + T(14)*c ,
                 T( 7)*s + T(15)*c , T( 8)*s + T(16)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RY0.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #endif

    vec4 = v4 ;
    RY0.apply( qclab::Op::Trans , 4 , vec4 ) ;
    check4 = {  T( 1)*c + T( 9)*s ,  T( 2)*c + T(10)*s ,
                T( 3)*c + T(11)*s ,  T( 4)*c + T(12)*s ,
                T( 5)*c + T(13)*s ,  T( 6)*c + T(14)*s ,
                T( 7)*c + T(15)*s ,  T( 8)*c + T(16)*s ,
               -T( 1)*s + T( 9)*c , -T( 2)*s + T(10)*c ,
               -T( 3)*s + T(11)*c , -T( 4)*s + T(12)*c ,
               -T( 5)*s + T(13)*c , -T( 6)*s + T(14)*c ,
               -T( 7)*s + T(15)*c , -T( 8)*s + T(16)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RY0.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #endif

    vec4 = v4 ;
    RY0.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RY0.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #endif

    vec4 = v4 ;
    RY1.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = { T( 1)*c - T( 5)*s , T( 2)*c - T( 6)*s ,
               T( 3)*c - T( 7)*s , T( 4)*c - T( 8)*s ,
               T( 1)*s + T( 5)*c , T( 2)*s + T( 6)*c ,
               T( 3)*s + T( 7)*c , T( 4)*s + T( 8)*c ,
               T( 9)*c - T(13)*s , T(10)*c - T(14)*s ,
               T(11)*c - T(15)*s , T(12)*c - T(16)*s ,
               T( 9)*s + T(13)*c , T(10)*s + T(14)*c ,
               T(11)*s + T(15)*c , T(12)*s + T(16)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RY1.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #endif

    vec4 = v4 ;
    RY1.apply( qclab::Op::Trans , 4 , vec4 ) ;
    check4 = {  T( 1)*c + T( 5)*s ,  T( 2)*c + T( 6)*s ,
                T( 3)*c + T( 7)*s ,  T( 4)*c + T( 8)*s ,
               -T( 1)*s + T( 5)*c , -T( 2)*s + T( 6)*c ,
               -T( 3)*s + T( 7)*c , -T( 4)*s + T( 8)*c ,
                T( 9)*c + T(13)*s ,  T(10)*c + T(14)*s ,
                T(11)*c + T(15)*s ,  T(12)*c + T(16)*s ,
               -T( 9)*s + T(13)*c , -T(10)*s + T(14)*c ,
               -T(11)*s + T(15)*c , -T(12)*s + T(16)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RY1.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #endif

    vec4 = v4 ;
    RY1.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RY1.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #endif

    vec4 = v4 ;
    RY2.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = { T( 1)*c - T( 3)*s , T( 2)*c - T( 4)*s ,
               T( 1)*s + T( 3)*c , T( 2)*s + T( 4)*c ,
               T( 5)*c - T( 7)*s , T( 6)*c - T( 8)*s ,
               T( 5)*s + T( 7)*c , T( 6)*s + T( 8)*c ,
               T( 9)*c - T(11)*s , T(10)*c - T(12)*s ,
               T( 9)*s + T(11)*c , T(10)*s + T(12)*c ,
               T(13)*c - T(15)*s , T(14)*c - T(16)*s ,
               T(13)*s + T(15)*c , T(14)*s + T(16)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RY2.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #endif

    vec4 = v4 ;
    RY2.apply( qclab::Op::Trans , 4 , vec4 ) ;
    check4 = {  T( 1)*c + T( 3)*s ,  T( 2)*c + T( 4)*s ,
               -T( 1)*s + T( 3)*c , -T( 2)*s + T( 4)*c ,
                T( 5)*c + T( 7)*s ,  T( 6)*c + T( 8)*s ,
               -T( 5)*s + T( 7)*c , -T( 6)*s + T( 8)*c ,
                T( 9)*c + T(11)*s ,  T(10)*c + T(12)*s ,
               -T( 9)*s + T(11)*c , -T(10)*s + T(12)*c ,
                T(13)*c + T(15)*s ,  T(14)*c + T(16)*s ,
               -T(13)*s + T(15)*c , -T(14)*s + T(16)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RY2.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #endif

    vec4 = v4 ;
    RY2.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RY2.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #endif

    vec4 = v4 ;
    RY3.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = { T( 1)*c - T( 2)*s , T( 1)*s + T( 2)*c ,
               T( 3)*c - T( 4)*s , T( 3)*s + T( 4)*c ,
               T( 5)*c - T( 6)*s , T( 5)*s + T( 6)*c ,
               T( 7)*c - T( 8)*s , T( 7)*s + T( 8)*c ,
               T( 9)*c - T(10)*s , T( 9)*s + T(10)*c ,
               T(11)*c - T(12)*s , T(11)*s + T(12)*c ,
               T(13)*c - T(14)*s , T(13)*s + T(14)*c ,
               T(15)*c - T(16)*s , T(15)*s + T(16)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RY3.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #endif

    vec4 = v4 ;
    RY3.apply( qclab::Op::Trans , 4 , vec4 ) ;
    check4 = { T( 1)*c + T( 2)*s , -T( 1)*s + T( 2)*c ,
               T( 3)*c + T( 4)*s , -T( 3)*s + T( 4)*c ,
               T( 5)*c + T( 6)*s , -T( 5)*s + T( 6)*c ,
               T( 7)*c + T( 8)*s , -T( 7)*s + T( 8)*c ,
               T( 9)*c + T(10)*s , -T( 9)*s + T(10)*c ,
               T(11)*c + T(12)*s , -T(11)*s + T(12)*c ,
               T(13)*c + T(14)*s , -T(13)*s + T(14)*c ,
               T(15)*c + T(16)*s , -T(15)*s + T(16)*c } ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RY3.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #endif

    vec4 = v4 ;
    RY3.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RY3.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( check_qclab_qgates_RotationY( vec4 , check4 ) ) ;
  #endif
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

