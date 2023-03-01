#include <gtest/gtest.h>
#include "qclab/qgates/RotationZ.hpp"

template <typename T>
void test_qclab_qgates_RotationZ() {

  using R = qclab::real_t< T > ;
  const R pi = 4 * std::atan(1) ;
  const R tol = 10 * std::numeric_limits< R >::epsilon() ;

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

    // update(rot)
    qclab::QRotation< R >  new_rot( 1 ) ;
    Rz.update( new_rot ) ;
    EXPECT_NEAR( Rz.theta() , 1 , tol ) ;
    EXPECT_NEAR( Rz.cos() , std::cos( 0.5 ) , tol ) ;
    EXPECT_NEAR( Rz.sin() , std::sin( 0.5 ) , tol ) ;

    // update(theta)
    Rz.update( pi/2 ) ;
    EXPECT_NEAR( Rz.theta() , pi/2 , tol ) ;
    EXPECT_NEAR( Rz.cos() , std::cos( pi/4 ) , tol ) ;
    EXPECT_NEAR( Rz.sin() , std::sin( pi/4 ) , tol ) ;

    // matrix
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    EXPECT_NEAR( std::real( Rz.matrix()(0,0) ) ,  cos , tol ) ;
    EXPECT_NEAR( std::imag( Rz.matrix()(0,0) ) , -sin , tol ) ;
    EXPECT_EQ( Rz.matrix()(1,0) , T(0.0) ) ;
    EXPECT_EQ( Rz.matrix()(0,1) , T(0.0) ) ;
    EXPECT_NEAR( std::real( Rz.matrix()(1,1) ) ,  cos , tol ) ;
    EXPECT_NEAR( std::imag( Rz.matrix()(1,1) ) ,  sin , tol ) ;

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
    EXPECT_NEAR( Rz.theta() , 2*(pi/3) , tol ) ;

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

    EXPECT_NEAR( Rz.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Rz.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Rz.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    qclab::qgates::RotationZ< T >  Rz( 5 , pi/2 ) ;

    EXPECT_EQ( Rz.nbQubits() , 1 ) ;                    // nbQubits
    EXPECT_FALSE( Rz.fixed() ) ;                        // fixed
    EXPECT_FALSE( Rz.controlled() ) ;                   // controlled
    EXPECT_EQ( Rz.qubit() , 5 ) ;                       // qubit

    EXPECT_NEAR( Rz.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Rz.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Rz.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    qclab::qgates::RotationZ< T >  Rz( 5 , pi/2 , true ) ;

    EXPECT_EQ( Rz.nbQubits() , 1 ) ;                    // nbQubits
    EXPECT_TRUE( Rz.fixed() ) ;                         // fixed
    EXPECT_FALSE( Rz.controlled() ) ;                   // controlled
    EXPECT_EQ( Rz.qubit() , 5 ) ;                       // qubit

    EXPECT_NEAR( Rz.theta() , pi/2 , tol ) ;            // theta
    EXPECT_NEAR( Rz.cos() , std::cos( pi/4 ) , tol ) ;  // cos
    EXPECT_NEAR( Rz.sin() , std::sin( pi/4 ) , tol ) ;  // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationZ< T >  Rz( cos , sin ) ;

    EXPECT_EQ( Rz.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_FALSE( Rz.fixed() ) ;              // fixed
    EXPECT_FALSE( Rz.controlled() ) ;         // controlled
    EXPECT_EQ( Rz.qubit() , 0 ) ;             // qubit

    EXPECT_NEAR( Rz.theta() , pi/2 , tol ) ;  // theta
    EXPECT_NEAR( Rz.cos() , cos , tol ) ;     // cos
    EXPECT_NEAR( Rz.sin() , sin , tol ) ;     // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationZ< T >  Rz( 5 , cos , sin ) ;

    EXPECT_EQ( Rz.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_FALSE( Rz.fixed() ) ;              // fixed
    EXPECT_FALSE( Rz.controlled() ) ;         // controlled
    EXPECT_EQ( Rz.qubit() , 5 ) ;             // qubit

    EXPECT_NEAR( Rz.theta() , pi/2 , tol ) ;  // theta
    EXPECT_NEAR( Rz.cos() , cos , tol ) ;     // cos
    EXPECT_NEAR( Rz.sin() , sin , tol ) ;     // sin
  }

  {
    const R cos = std::cos( pi/4 ) ;
    const R sin = std::sin( pi/4 ) ;
    qclab::qgates::RotationZ< T >  Rz( 5 , cos , sin , true ) ;

    EXPECT_EQ( Rz.nbQubits() , 1 ) ;          // nbQubits
    EXPECT_TRUE( Rz.fixed() ) ;               // fixed
    EXPECT_FALSE( Rz.controlled() ) ;         // controlled
    EXPECT_EQ( Rz.qubit() , 5 ) ;             // qubit

    EXPECT_NEAR( Rz.theta() , pi/2 , tol ) ;  // theta
    EXPECT_NEAR( Rz.cos() , cos , tol ) ;     // cos
    EXPECT_NEAR( Rz.sin() , sin , tol ) ;     // sin
  }

  {
    const R theta1 = pi/2 ;
    const qclab::QAngle< R >       angle1( theta1 ) ;
    const qclab::QRotation< R >    rot1( angle1 ) ;
    qclab::qgates::RotationZ< T >  R1( rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >       angle2( theta2 ) ;
    const qclab::QRotation< R >    rot2( angle2 ) ;
    qclab::qgates::RotationZ< T >  R2( rot2 ) ;

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
    qclab::qgates::RotationZ< T >  R1( rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >       angle2( theta2 ) ;
    const qclab::QRotation< R >    rot2( angle2 ) ;
    qclab::qgates::RotationZ< T >  R2( rot2 ) ;

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
    qclab::qgates::RotationZ< T >  R1( rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >       angle2( theta2 ) ;
    const qclab::QRotation< R >    rot2( angle2 ) ;
    qclab::qgates::RotationZ< T >  R2( rot2 ) ;

    // operator *
    qclab::qgates::RotationZ< T > R12 = R1 * R2 ;
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
    qclab::qgates::RotationZ< T >  R1( rot1 ) ;

    const R theta2 = pi/3 ;
    const qclab::QAngle< R >       angle2( theta2 ) ;
    const qclab::QRotation< R >    rot2( angle2 ) ;
    qclab::qgates::RotationZ< T >  R2( rot2 ) ;

    // operator /
    qclab::qgates::RotationZ< T > R12 = R1 / R2 ;
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
    qclab::qgates::RotationZ< T >  R1( rot1 ) ;
    R cos   = R1.cos() ;
    R sin   = R1.sin() ;
    R theta = R1.theta() ;

    // inv
    qclab::qgates::RotationZ< T >  R2 = R1.inv() ;
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

    const T c1 = T( std::cos( pi/8 ) , -std::sin( pi/8 ) ) ;
    const T c2 = T( std::cos( pi/8 ) ,  std::sin( pi/8 ) ) ;

    qclab::qgates::RotationZ< T >  RZ0( 0 , pi/4 ) ;
    qclab::qgates::RotationZ< T >  RZ1( 1 , pi/4 ) ;
    qclab::qgates::RotationZ< T >  RZ2( 2 , pi/4 ) ;
    qclab::qgates::RotationZ< T >  RZ3( 3 , pi/4 ) ;

    // nbQubits = 1
    auto vec1 = v1 ;
    RZ0.apply( qclab::Op::NoTrans , 1 , vec1 ) ;
    V check1 = { T(1)*c1 , T(2)*c2 } ;
    EXPECT_TRUE( vec1 == check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ; T* vec1_ = vec1.data() ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { RZ0.apply_device( qclab::Op::NoTrans , 1 , vec1_ ) ; }
    EXPECT_TRUE( vec1 == check1 ) ;
  #endif

    vec1 = v1 ;
    RZ0.apply( qclab::Op::Trans , 1 , vec1 ) ;
    EXPECT_TRUE( vec1 == check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { RZ0.apply_device( qclab::Op::Trans , 1 , vec1_ ) ; }
    EXPECT_TRUE( vec1 == check1 ) ;
  #endif

    vec1 = v1 ;
    RZ0.apply( qclab::Op::ConjTrans , 1 , vec1 ) ;
    check1 = { T(1)*c2 , T(2)*c1 } ;
    EXPECT_TRUE( vec1 == check1 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec1 = v1 ;
    #pragma omp target data map(tofrom:vec1_[0:2])
    { RZ0.apply_device( qclab::Op::ConjTrans , 1 , vec1_ ) ; }
    EXPECT_TRUE( vec1 == check1 ) ;
  #endif

    // nbQubits = 2
    auto vec2 = v2 ;
    RZ0.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { T(1)*c1 , T(2)*c1 , T(3)*c2 , T(4)*c2 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { RZ0.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    RZ0.apply( qclab::Op::Trans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { RZ0.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    RZ0.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check2 = { T(1)*c2 , T(2)*c2 , T(3)*c1 , T(4)*c1 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { RZ0.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    RZ1.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    check2 = { T(1)*c1 , T(2)*c2 , T(3)*c1 , T(4)*c2 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { RZ1.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    RZ1.apply( qclab::Op::Trans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { RZ1.apply_device( qclab::Op::Trans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    RZ1.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    check2 = { T(1)*c2 , T(2)*c1 , T(3)*c2 , T(4)*c1 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { RZ1.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    RZ0.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { T(1)*c1 , T(2)*c1 , T(3)*c1 , T(4)*c1 ,
                 T(5)*c2 , T(6)*c2 , T(7)*c2 , T(8)*c2 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RZ0.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    RZ0.apply( qclab::Op::Trans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RZ0.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    RZ0.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { T(1)*c2 , T(2)*c2 , T(3)*c2 , T(4)*c2 ,
               T(5)*c1 , T(6)*c1 , T(7)*c1 , T(8)*c1 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RZ0.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    RZ1.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T(1)*c1 , T(2)*c1 , T(3)*c2 , T(4)*c2 ,
               T(5)*c1 , T(6)*c1 , T(7)*c2 , T(8)*c2 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RZ1.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    RZ1.apply( qclab::Op::Trans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RZ1.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    RZ1.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { T(1)*c2 , T(2)*c2 , T(3)*c1 , T(4)*c1 ,
               T(5)*c2 , T(6)*c2 , T(7)*c1 , T(8)*c1 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RZ1.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    RZ2.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { T(1)*c1 , T(2)*c2 , T(3)*c1 , T(4)*c2 ,
               T(5)*c1 , T(6)*c2 , T(7)*c1 , T(8)*c2 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RZ2.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    RZ2.apply( qclab::Op::Trans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RZ2.apply_device( qclab::Op::Trans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    RZ2.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    check3 = { T(1)*c2 , T(2)*c1 , T(3)*c2 , T(4)*c1 ,
               T(5)*c2 , T(6)*c1 , T(7)*c2 , T(8)*c1 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { RZ2.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    RZ0.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = { T( 1)*c1 , T( 2)*c1 , T( 3)*c1 , T( 4)*c1 ,
                 T( 5)*c1 , T( 6)*c1 , T( 7)*c1 , T( 8)*c1 ,
                 T( 9)*c2 , T(10)*c2 , T(11)*c2 , T(12)*c2 ,
                 T(13)*c2 , T(14)*c2 , T(15)*c2 , T(16)*c2 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RZ0.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    RZ0.apply( qclab::Op::Trans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RZ0.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    RZ0.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { T( 1)*c2 , T( 2)*c2 , T( 3)*c2 , T( 4)*c2 ,
               T( 5)*c2 , T( 6)*c2 , T( 7)*c2 , T( 8)*c2 ,
               T( 9)*c1 , T(10)*c1 , T(11)*c1 , T(12)*c1 ,
               T(13)*c1 , T(14)*c1 , T(15)*c1 , T(16)*c1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RZ0.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    RZ1.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = { T( 1)*c1 , T( 2)*c1 , T( 3)*c1 , T( 4)*c1 ,
               T( 5)*c2 , T( 6)*c2 , T( 7)*c2 , T( 8)*c2 ,
               T( 9)*c1 , T(10)*c1 , T(11)*c1 , T(12)*c1 ,
               T(13)*c2 , T(14)*c2 , T(15)*c2 , T(16)*c2 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RZ1.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    RZ1.apply( qclab::Op::Trans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RZ1.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    RZ1.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { T( 1)*c2 , T( 2)*c2 , T( 3)*c2 , T( 4)*c2 ,
               T( 5)*c1 , T( 6)*c1 , T( 7)*c1 , T( 8)*c1 ,
               T( 9)*c2 , T(10)*c2 , T(11)*c2 , T(12)*c2 ,
               T(13)*c1 , T(14)*c1 , T(15)*c1 , T(16)*c1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RZ1.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    RZ2.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = { T( 1)*c1 , T( 2)*c1 , T( 3)*c2 , T( 4)*c2 ,
               T( 5)*c1 , T( 6)*c1 , T( 7)*c2 , T( 8)*c2 ,
               T( 9)*c1 , T(10)*c1 , T(11)*c2 , T(12)*c2 ,
               T(13)*c1 , T(14)*c1 , T(15)*c2 , T(16)*c2 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RZ2.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    RZ2.apply( qclab::Op::Trans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RZ2.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    RZ2.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { T( 1)*c2 , T( 2)*c2 , T( 3)*c1 , T( 4)*c1 ,
               T( 5)*c2 , T( 6)*c2 , T( 7)*c1 , T( 8)*c1 ,
               T( 9)*c2 , T(10)*c2 , T(11)*c1 , T(12)*c1 ,
               T(13)*c2 , T(14)*c2 , T(15)*c1 , T(16)*c1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RZ2.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    RZ3.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    check4 = { T( 1)*c1 , T( 2)*c2 , T( 3)*c1 , T( 4)*c2 ,
               T( 5)*c1 , T( 6)*c2 , T( 7)*c1 , T( 8)*c2 ,
               T( 9)*c1 , T(10)*c2 , T(11)*c1 , T(12)*c2 ,
               T(13)*c1 , T(14)*c2 , T(15)*c1 , T(16)*c2 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RZ3.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    RZ3.apply( qclab::Op::Trans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RZ3.apply_device( qclab::Op::Trans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    RZ3.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    check4 = { T( 1)*c2 , T( 2)*c1 , T( 3)*c2 , T( 4)*c1 ,
               T( 5)*c2 , T( 6)*c1 , T( 7)*c2 , T( 8)*c1 ,
               T( 9)*c2 , T(10)*c1 , T(11)*c2 , T(12)*c1 ,
               T(13)*c2 , T(14)*c1 , T(15)*c2 , T(16)*c1 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { RZ3.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif
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

