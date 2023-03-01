#include <gtest/gtest.h>
#include "qclab/qgates/CZ.hpp"
#include <numeric>

template <typename T>
void test_qclab_qgates_CZ() {

  //
  // CZs |1> controlled
  //
  {
    qclab::qgates::CZ< T >  cz ;

    EXPECT_EQ( cz.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cz.fixed() ) ;           // fixed
    EXPECT_TRUE( cz.controlled() ) ;      // controlled
    EXPECT_EQ( cz.control() , 0 ) ;       // control
    EXPECT_EQ( cz.target() , 1 ) ;        // target
    EXPECT_EQ( cz.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CZ_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 0 , 0 ,-1 ) ;
    EXPECT_TRUE( cz.matrix() == CZ_check ) ;

    // qubit
    EXPECT_EQ( cz.qubit() , 0 ) ;

    // qubits
    auto qubits = cz.qubits() ;
    EXPECT_EQ( qubits.size() , 2 ) ;
    EXPECT_EQ( qubits[0] , 0 ) ;
    EXPECT_EQ( qubits[1] , 1 ) ;
    int qnew[] = { 5 , 3 } ;
    cz.setQubits( &qnew[0] ) ;
    EXPECT_EQ( cz.qubits()[0] , 3 ) ;
    EXPECT_EQ( cz.qubits()[1] , 5 ) ;
    qnew[0] = 0 ;
    qnew[1] = 1 ;
    cz.setQubits( &qnew[0] ) ;

    // print
    cz.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cz.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cz q[0], q[1];\n" ) ;
    std::cout << qasm.str() ;

    // gate
    qclab::qgates::PauliZ< T >  Z ;
    EXPECT_TRUE( *cz.gate() == Z ) ;
    EXPECT_TRUE( cz.gate()->matrix() == Z.matrix() ) ;

    // operators == and !=
    EXPECT_TRUE(  cz != Z ) ;
    EXPECT_FALSE( cz == Z ) ;
    qclab::qgates::CZ< T >  cz2 ;
    EXPECT_TRUE(  cz == cz2 ) ;
    EXPECT_FALSE( cz != cz2 ) ;

    // setControl, setTarget, setControlState
    cz.setControl( 3 ) ;
    EXPECT_EQ( cz.control() , 3 ) ;
    cz.setTarget( 5 ) ;
    EXPECT_EQ( cz.target() , 5 ) ;
    EXPECT_TRUE(  cz == cz2 ) ;
    EXPECT_FALSE( cz != cz2 ) ;

    cz.setControl( 4 ) ;
    EXPECT_EQ( cz.control() , 4 ) ;
    cz.setTarget( 1 ) ;
    EXPECT_EQ( cz.target() , 1 ) ;
    EXPECT_TRUE(  cz != cz2 ) ;
    EXPECT_FALSE( cz == cz2 ) ;

    qubits[0] = 1 ;
    qubits[1] = 2 ;
    cz.setQubits( &qubits[0] ) ;
    EXPECT_EQ( cz.control() , 1 ) ;
    EXPECT_EQ( cz.target() , 2 ) ;
    EXPECT_TRUE(  cz == cz2 ) ;
    EXPECT_FALSE( cz != cz2 ) ;

    cz.setControl( 0 ) ;
    EXPECT_EQ( cz.control() , 0 ) ;
    cz.setTarget( 1 ) ;
    EXPECT_EQ( cz.target() , 1 ) ;
    cz.setControlState( 0 ) ;
    EXPECT_EQ( cz.controlState() , 0 ) ;
    EXPECT_TRUE(  cz != cz2 ) ;
    EXPECT_FALSE( cz == cz2 ) ;
  }

  {
    qclab::qgates::CZ< T >  cz01( 0 , 1 ) ;

    EXPECT_EQ( cz01.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cz01.fixed() ) ;           // fixed
    EXPECT_TRUE( cz01.controlled() ) ;      // controlled
    EXPECT_EQ( cz01.control() , 0 ) ;       // control
    EXPECT_EQ( cz01.target() , 1 ) ;        // target
    EXPECT_EQ( cz01.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CZ_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 0 , 0 ,-1 ) ;
    EXPECT_TRUE( cz01.matrix() == CZ_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cz01.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cz q[0], q[1];\n" ) ;
  }

  {
    qclab::qgates::CZ< T >  cz35( 3 , 5 ) ;

    EXPECT_EQ( cz35.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cz35.fixed() ) ;           // fixed
    EXPECT_TRUE( cz35.controlled() ) ;      // controlled
    EXPECT_EQ( cz35.control() , 3 ) ;       // control
    EXPECT_EQ( cz35.target() , 5 ) ;        // target
    EXPECT_EQ( cz35.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CZ_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 0 , 0 ,-1 ) ;
    EXPECT_TRUE( cz35.matrix() == CZ_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cz35.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cz q[3], q[5];\n" ) ;
  }

  {
    qclab::qgates::CZ< T >  cz10( 1 , 0 ) ;

    EXPECT_EQ( cz10.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cz10.fixed() ) ;           // fixed
    EXPECT_TRUE( cz10.controlled() ) ;      // controlled
    EXPECT_EQ( cz10.control() , 1 ) ;       // control
    EXPECT_EQ( cz10.target() , 0 ) ;        // target
    EXPECT_EQ( cz10.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CZ_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 0 , 0 ,-1 ) ;
    EXPECT_TRUE( cz10.matrix() == CZ_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cz10.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cz q[1], q[0];\n" ) ;
  }

  {
    qclab::qgates::CZ< T >  cz53( 5 , 3 ) ;

    EXPECT_EQ( cz53.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cz53.fixed() ) ;           // fixed
    EXPECT_TRUE( cz53.controlled() ) ;      // controlled
    EXPECT_EQ( cz53.control() , 5 ) ;       // control
    EXPECT_EQ( cz53.target() , 3 ) ;        // target
    EXPECT_EQ( cz53.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CZ_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 0 , 0 ,-1 ) ;
    EXPECT_TRUE( cz53.matrix() == CZ_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cz53.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cz q[5], q[3];\n" ) ;
  }


  //
  // CZs |0> controlled
  //
  {
    qclab::qgates::CZ< T >  cz01( 0 , 1 , 0 ) ;

    EXPECT_EQ( cz01.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cz01.fixed() ) ;           // fixed
    EXPECT_TRUE( cz01.controlled() ) ;      // controlled
    EXPECT_EQ( cz01.control() , 0 ) ;       // control
    EXPECT_EQ( cz01.target() , 1 ) ;        // target
    EXPECT_EQ( cz01.controlState() , 0 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CZ_check( 1 , 0 , 0 , 0 ,
                                               0 ,-1 , 0 , 0 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 0 , 0 , 1 ) ;
    EXPECT_TRUE( cz01.matrix() == CZ_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cz01.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[0];\ncz q[0], q[1];\nx q[0];\n" ) ;
  }

  {
    qclab::qgates::CZ< T >  cz10( 1 , 0 , 0 ) ;

    EXPECT_EQ( cz10.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cz10.fixed() ) ;           // fixed
    EXPECT_TRUE( cz10.controlled() ) ;      // controlled
    EXPECT_EQ( cz10.control() , 1 ) ;       // control
    EXPECT_EQ( cz10.target() , 0 ) ;        // target
    EXPECT_EQ( cz10.controlState() , 0 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CZ_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 ,-1 , 0 ,
                                               0 , 0 , 0 , 1 ) ;
    EXPECT_TRUE( cz10.matrix() == CZ_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cz10.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[1];\ncz q[1], q[0];\nx q[1];\n" ) ;
  }

  using V = std::vector< T > ;

  const V v2 = { 1 , 2 , 3 , 4 } ;
  const V v3 = { 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 } ;
  V vt( 16 ) ; std::iota( vt.begin() , vt.end() , 1 ) ;
  const V v4 = vt ;
  vt = V( 32 ) ; std::iota( vt.begin() , vt.end() , 1 ) ;
  const V v5 = vt ;
  vt = V( 256 ) ; std::iota( vt.begin() , vt.end() , 1 ) ;
  const V v8 = vt ;

  //
  // CNOTs |1> controlled
  //

  // control < target
  {
    qclab::qgates::CZ< T >  cz01( 0 , 1 ) ;
    qclab::qgates::CZ< T >  cz12( 1 , 2 ) ;
    qclab::qgates::CZ< T >  cz02( 0 , 2 ) ;
    qclab::qgates::CZ< T >  cz13( 1 , 3 ) ;
    qclab::qgates::CZ< T >  cz26( 2 , 6 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cz01.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 1 , 2 , 3 , -4 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cz01.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cz01.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cz01.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    cz01.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 1 , 2 , 3 , 4 , 5 , 6 , -7 , -8 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz01.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz01.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz01.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz12.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 1 , 2 , 3 , -4 , 5 , 6 , 7 , -8 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz12.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz12.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz12.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz02.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 1 , 2 , 3 , 4 , 5 , -6 , 7 , -8 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz02.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz02.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz02.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    cz12.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = {   1 ,   2 ,   3 ,   4 ,   5 ,   6 ,  -7 ,  -8 ,
                   9 ,  10 ,  11 ,  12 ,  13 ,  14 , -15 , -16 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cz12.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cz12.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cz12.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    cz13.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = {   1 ,   2 ,   3 ,   4 ,   5 ,   6 ,   7 ,   8 ,
                   9 ,  10 , -11 , -12 ,  13 ,  14 , -15 , -16 ,
                  17 ,  18 ,  19 ,  20 ,  21 ,  22 ,  23 ,  24 ,
                  25 ,  26 , -27 , -28 ,  29 ,  30 , -31 , -32 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cz13.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cz13.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cz13.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    cz26.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    V check8 = {    1 ,    2 ,    3 ,    4 ,    5 ,    6 ,    7 ,    8 ,
                    9 ,   10 ,   11 ,   12 ,   13 ,   14 ,   15 ,   16 ,
                   17 ,   18 ,   19 ,   20 ,   21 ,   22 ,   23 ,   24 ,
                   25 ,   26 ,   27 ,   28 ,   29 ,   30 ,   31 ,   32 ,
                   33 ,   34 ,  -35 ,  -36 ,   37 ,   38 ,  -39 ,  -40 ,
                   41 ,   42 ,  -43 ,  -44 ,   45 ,   46 ,  -47 ,  -48 ,
                   49 ,   50 ,  -51 ,  -52 ,   53 ,   54 ,  -55 ,  -56 ,
                   57 ,   58 ,  -59 ,  -60 ,   61 ,   62 ,  -63 ,  -64 ,
                   65 ,   66 ,   67 ,   68 ,   69 ,   70 ,   71 ,   72 ,
                   73 ,   74 ,   75 ,   76 ,   77 ,   78 ,   79 ,   80 ,
                   81 ,   82 ,   83 ,   84 ,   85 ,   86 ,   87 ,   88 ,
                   89 ,   90 ,   91 ,   92 ,   93 ,   94 ,   95 ,   96 ,
                   97 ,   98 ,  -99 , -100 ,  101 ,  102 , -103 , -104 ,
                  105 ,  106 , -107 , -108 ,  109 ,  110 , -111 , -112 ,
                  113 ,  114 , -115 , -116 ,  117 ,  118 , -119 , -120 ,
                  121 ,  122 , -123 , -124 ,  125 ,  126 , -127 , -128 ,
                  129 ,  130 ,  131 ,  132 ,  133 ,  134 ,  135 ,  136 ,
                  137 ,  138 ,  139 ,  140 ,  141 ,  142 ,  143 ,  144 ,
                  145 ,  146 ,  147 ,  148 ,  149 ,  150 ,  151 ,  152 ,
                  153 ,  154 ,  155 ,  156 ,  157 ,  158 ,  159 ,  160 ,
                  161 ,  162 , -163 , -164 ,  165 ,  166 , -167 , -168 ,
                  169 ,  170 , -171 , -172 ,  173 ,  174 , -175 , -176 ,
                  177 ,  178 , -179 , -180 ,  181 ,  182 , -183 , -184 ,
                  185 ,  186 , -187 , -188 ,  189 ,  190 , -191 , -192 ,
                  193 ,  194 ,  195 ,  196 ,  197 ,  198 ,  199 ,  200 ,
                  201 ,  202 ,  203 ,  204 ,  205 ,  206 ,  207 ,  208 ,
                  209 ,  210 ,  211 ,  212 ,  213 ,  214 ,  215 ,  216 ,
                  217 ,  218 ,  219 ,  220 ,  221 ,  222 ,  223 ,  224 ,
                  225 ,  226 , -227 , -228 ,  229 ,  230 , -231 , -232 ,
                  233 ,  234 , -235 , -236 ,  237 ,  238 , -239 , -240 ,
                  241 ,  242 , -243 , -244 ,  245 ,  246 , -247 , -248 ,
                  249 ,  250 , -251 , -252 ,  253 ,  254 , -255 , -256 } ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cz26.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cz26.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cz26.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

  // control > target
  {
    qclab::qgates::CZ< T >  cz10( 1 , 0 ) ;
    qclab::qgates::CZ< T >  cz21( 2 , 1 ) ;
    qclab::qgates::CZ< T >  cz20( 2 , 0 ) ;
    qclab::qgates::CZ< T >  cz31( 3 , 1 ) ;
    qclab::qgates::CZ< T >  cz62( 6 , 2 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cz10.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 1 , 2 , 3 , -4 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cz10.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cz10.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cz10.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    cz10.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 1 , 2 , 3 , 4 , 5 , 6 , -7 , -8 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz10.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz10.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz10.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz21.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 1 , 2 , 3 , -4 , 5 , 6 , 7 , -8 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz21.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz21.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz21.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz20.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 1 , 2 , 3 , 4 , 5 , -6 , 7 , -8 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz20.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz20.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz20.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    cz21.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = {   1 ,   2 ,   3 ,   4 ,   5 ,   6 ,  -7 ,  -8 ,
                   9 ,  10 ,  11 ,  12 ,  13 ,  14 , -15 , -16 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cz21.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cz21.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cz21.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    cz31.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = {   1 ,   2 ,   3 ,   4 ,   5 ,   6 ,   7 ,   8 ,
                   9 ,  10 , -11 , -12 ,  13 ,  14 , -15 , -16 ,
                  17 ,  18 ,  19 ,  20 ,  21 ,  22 ,  23 ,  24 ,
                  25 ,  26 , -27 , -28 ,  29 ,  30 , -31 , -32 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cz31.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cz31.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cz31.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    cz62.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    V check8 = {    1 ,    2 ,    3 ,    4 ,    5 ,    6 ,    7 ,    8 ,
                    9 ,   10 ,   11 ,   12 ,   13 ,   14 ,   15 ,   16 ,
                   17 ,   18 ,   19 ,   20 ,   21 ,   22 ,   23 ,   24 ,
                   25 ,   26 ,   27 ,   28 ,   29 ,   30 ,   31 ,   32 ,
                   33 ,   34 ,  -35 ,  -36 ,   37 ,   38 ,  -39 ,  -40 ,
                   41 ,   42 ,  -43 ,  -44 ,   45 ,   46 ,  -47 ,  -48 ,
                   49 ,   50 ,  -51 ,  -52 ,   53 ,   54 ,  -55 ,  -56 ,
                   57 ,   58 ,  -59 ,  -60 ,   61 ,   62 ,  -63 ,  -64 ,
                   65 ,   66 ,   67 ,   68 ,   69 ,   70 ,   71 ,   72 ,
                   73 ,   74 ,   75 ,   76 ,   77 ,   78 ,   79 ,   80 ,
                   81 ,   82 ,   83 ,   84 ,   85 ,   86 ,   87 ,   88 ,
                   89 ,   90 ,   91 ,   92 ,   93 ,   94 ,   95 ,   96 ,
                   97 ,   98 ,  -99 , -100 ,  101 ,  102 , -103 , -104 ,
                  105 ,  106 , -107 , -108 ,  109 ,  110 , -111 , -112 ,
                  113 ,  114 , -115 , -116 ,  117 ,  118 , -119 , -120 ,
                  121 ,  122 , -123 , -124 ,  125 ,  126 , -127 , -128 ,
                  129 ,  130 ,  131 ,  132 ,  133 ,  134 ,  135 ,  136 ,
                  137 ,  138 ,  139 ,  140 ,  141 ,  142 ,  143 ,  144 ,
                  145 ,  146 ,  147 ,  148 ,  149 ,  150 ,  151 ,  152 ,
                  153 ,  154 ,  155 ,  156 ,  157 ,  158 ,  159 ,  160 ,
                  161 ,  162 , -163 , -164 ,  165 ,  166 , -167 , -168 ,
                  169 ,  170 , -171 , -172 ,  173 ,  174 , -175 , -176 ,
                  177 ,  178 , -179 , -180 ,  181 ,  182 , -183 , -184 ,
                  185 ,  186 , -187 , -188 ,  189 ,  190 , -191 , -192 ,
                  193 ,  194 ,  195 ,  196 ,  197 ,  198 ,  199 ,  200 ,
                  201 ,  202 ,  203 ,  204 ,  205 ,  206 ,  207 ,  208 ,
                  209 ,  210 ,  211 ,  212 ,  213 ,  214 ,  215 ,  216 ,
                  217 ,  218 ,  219 ,  220 ,  221 ,  222 ,  223 ,  224 ,
                  225 ,  226 , -227 , -228 ,  229 ,  230 , -231 , -232 ,
                  233 ,  234 , -235 , -236 ,  237 ,  238 , -239 , -240 ,
                  241 ,  242 , -243 , -244 ,  245 ,  246 , -247 , -248 ,
                  249 ,  250 , -251 , -252 ,  253 ,  254 , -255 , -256 } ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cz62.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cz62.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cz62.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

  //
  // CNOTs |0> controlled
  //

  // control < target
  {
    qclab::qgates::CZ< T >  cz01( 0 , 1 , 0 ) ;
    qclab::qgates::CZ< T >  cz12( 1 , 2 , 0 ) ;
    qclab::qgates::CZ< T >  cz02( 0 , 2 , 0 ) ;
    qclab::qgates::CZ< T >  cz13( 1 , 3 , 0 ) ;
    qclab::qgates::CZ< T >  cz26( 2 , 6 , 0 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cz01.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 1 , -2 , 3 , 4 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cz01.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cz01.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cz01.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    cz01.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 1 , 2 , -3 , -4 , 5 , 6 , 7 , 8 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz01.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz01.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz01.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz12.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 1 , -2 , 3 , 4 , 5 , -6 , 7 , 8 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz12.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz12.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz12.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz02.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 1 , -2 , 3 , -4 , 5 , 6 , 7 , 8 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz02.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz02.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz02.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    cz12.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = {   1 ,   2 ,  -3 ,  -4 ,   5 ,   6 ,   7 ,   8 ,
                   9 ,  10 , -11 , -12 ,  13 ,  14 ,  15 ,  16 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cz12.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cz12.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cz12.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    cz13.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = {   1 ,   2 ,  -3 ,  -4 ,   5 ,   6 ,  -7 ,  -8 ,
                   9 ,  10 ,  11 ,  12 ,  13 ,  14 ,  15 ,  16 ,
                  17 ,  18 , -19 , -20 ,  21 ,  22 , -23 , -24 ,
                  25 ,  26 ,  27 ,  28 ,  29 ,  30 ,  31 ,  32 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cz13.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cz13.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cz13.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    cz26.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    V check8 = {    1 ,    2 ,   -3 ,   -4 ,    5 ,    6 ,   -7 ,   -8 ,
                    9 ,   10 ,  -11 ,  -12 ,   13 ,   14 ,  -15 ,  -16 ,
                   17 ,   18 ,  -19 ,  -20 ,   21 ,   22 ,  -23 ,  -24 ,
                   25 ,   26 ,  -27 ,  -28 ,   29 ,   30 ,  -31 ,  -32 ,
                   33 ,   34 ,   35 ,   36 ,   37 ,   38 ,   39 ,   40 ,
                   41 ,   42 ,   43 ,   44 ,   45 ,   46 ,   47 ,   48 ,
                   49 ,   50 ,   51 ,   52 ,   53 ,   54 ,   55 ,   56 ,
                   57 ,   58 ,   59 ,   60 ,   61 ,   62 ,   63 ,   64 ,
                   65 ,   66 ,  -67 ,  -68 ,   69 ,   70 ,  -71 ,  -72 ,
                   73 ,   74 ,  -75 ,  -76 ,   77 ,   78 ,  -79 ,  -80 ,
                   81 ,   82 ,  -83 ,  -84 ,   85 ,   86 ,  -87 ,  -88 ,
                   89 ,   90 ,  -91 ,  -92 ,   93 ,   94 ,  -95 ,  -96 ,
                   97 ,   98 ,   99 ,  100 ,  101 ,  102 ,  103 ,  104 ,
                  105 ,  106 ,  107 ,  108 ,  109 ,  110 ,  111 ,  112 ,
                  113 ,  114 ,  115 ,  116 ,  117 ,  118 ,  119 ,  120 ,
                  121 ,  122 ,  123 ,  124 ,  125 ,  126 ,  127 ,  128 ,
                  129 ,  130 , -131 , -132 ,  133 ,  134 , -135 , -136 ,
                  137 ,  138 , -139 , -140 ,  141 ,  142 , -143 , -144 ,
                  145 ,  146 , -147 , -148 ,  149 ,  150 , -151 , -152 ,
                  153 ,  154 , -155 , -156 ,  157 ,  158 , -159 , -160 ,
                  161 ,  162 ,  163 ,  164 ,  165 ,  166 ,  167 ,  168 ,
                  169 ,  170 ,  171 ,  172 ,  173 ,  174 ,  175 ,  176 ,
                  177 ,  178 ,  179 ,  180 ,  181 ,  182 ,  183 ,  184 ,
                  185 ,  186 ,  187 ,  188 ,  189 ,  190 ,  191 ,  192 ,
                  193 ,  194 , -195 , -196 ,  197 ,  198 , -199 , -200 ,
                  201 ,  202 , -203 , -204 ,  205 ,  206 , -207 , -208 ,
                  209 ,  210 , -211 , -212 ,  213 ,  214 , -215 , -216 ,
                  217 ,  218 , -219 , -220 ,  221 ,  222 , -223 , -224 ,
                  225 ,  226 ,  227 ,  228 ,  229 ,  230 ,  231 ,  232 ,
                  233 ,  234 ,  235 ,  236 ,  237 ,  238 ,  239 ,  240 ,
                  241 ,  242 ,  243 ,  244 ,  245 ,  246 ,  247 ,  248 ,
                  249 ,  250 ,  251 ,  252 ,  253 ,  254 ,  255 ,  256 } ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cz26.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cz26.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cz26.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

  // control > target
  {
    qclab::qgates::CZ< T >  cz10( 1 , 0 , 0 ) ;
    qclab::qgates::CZ< T >  cz21( 2 , 1 , 0 ) ;
    qclab::qgates::CZ< T >  cz20( 2 , 0 , 0 ) ;
    qclab::qgates::CZ< T >  cz31( 3 , 1 , 0 ) ;
    qclab::qgates::CZ< T >  cz62( 6 , 2 , 0 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cz10.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 1 , 2 , -3 , 4 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cz10.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cz10.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cz10.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    cz10.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 1 , 2 , 3 , 4 , -5 , -6 , 7 , 8 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz10.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz10.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz10.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz21.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 1 , 2 , -3 , 4 , 5 , 6 , -7 , 8 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz21.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz21.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz21.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz20.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 1 , 2 , 3 , 4 , -5 , 6 , -7 , 8 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz20.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cz20.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cz20.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    cz21.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = {   1 ,   2 ,   3 ,   4 ,  -5 ,  -6 ,   7 ,   8 ,
                   9 ,  10 ,  11 ,  12 , -13 , -14 ,  15 ,  16 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cz21.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cz21.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cz21.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    cz31.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = {   1 ,   2 ,   3 ,   4 ,   5 ,   6 ,   7 ,   8 ,
                  -9 , -10 ,  11 ,  12 , -13 , -14 ,  15 ,  16 ,
                  17 ,  18 ,  19 ,  20 ,  21 ,  22 ,  23 ,  24 ,
                 -25 , -26 ,  27 ,  28 , -29 , -30 ,  31 ,  32 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cz31.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cz31.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cz31.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    cz62.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    V check8 = {    1 ,    2 ,    3 ,    4 ,    5 ,    6 ,    7 ,    8 ,
                    9 ,   10 ,   11 ,   12 ,   13 ,   14 ,   15 ,   16 ,
                   17 ,   18 ,   19 ,   20 ,   21 ,   22 ,   23 ,   24 ,
                   25 ,   26 ,   27 ,   28 ,   29 ,   30 ,   31 ,   32 ,
                  -33 ,  -34 ,   35 ,   36 ,  -37 ,  -38 ,   39 ,   40 ,
                  -41 ,  -42 ,   43 ,   44 ,  -45 ,  -46 ,   47 ,   48 ,
                  -49 ,  -50 ,   51 ,   52 ,  -53 ,  -54 ,   55 ,   56 ,
                  -57 ,  -58 ,   59 ,   60 ,  -61 ,  -62 ,   63 ,   64 ,
                   65 ,   66 ,   67 ,   68 ,   69 ,   70 ,   71 ,   72 ,
                   73 ,   74 ,   75 ,   76 ,   77 ,   78 ,   79 ,   80 ,
                   81 ,   82 ,   83 ,   84 ,   85 ,   86 ,   87 ,   88 ,
                   89 ,   90 ,   91 ,   92 ,   93 ,   94 ,   95 ,   96 ,
                  -97 ,  -98 ,   99 ,  100 , -101 , -102 ,  103 ,  104 ,
                 -105 , -106 ,  107 ,  108 , -109 , -110 ,  111 ,  112 ,
                 -113 , -114 ,  115 ,  116 , -117 , -118 ,  119 ,  120 ,
                 -121 , -122 ,  123 ,  124 , -125 , -126 ,  127 ,  128 ,
                  129 ,  130 ,  131 ,  132 ,  133 ,  134 ,  135 ,  136 ,
                  137 ,  138 ,  139 ,  140 ,  141 ,  142 ,  143 ,  144 ,
                  145 ,  146 ,  147 ,  148 ,  149 ,  150 ,  151 ,  152 ,
                  153 ,  154 ,  155 ,  156 ,  157 ,  158 ,  159 ,  160 ,
                 -161 , -162 ,  163 ,  164 , -165 , -166 ,  167 ,  168 ,
                 -169 , -170 ,  171 ,  172 , -173 , -174 ,  175 ,  176 ,
                 -177 , -178 ,  179 ,  180 , -181 , -182 ,  183 ,  184 ,
                 -185 , -186 ,  187 ,  188 , -189 , -190 ,  191 ,  192 ,
                  193 ,  194 ,  195 ,  196 ,  197 ,  198 ,  199 ,  200 ,
                  201 ,  202 ,  203 ,  204 ,  205 ,  206 ,  207 ,  208 ,
                  209 ,  210 ,  211 ,  212 ,  213 ,  214 ,  215 ,  216 ,
                  217 ,  218 ,  219 ,  220 ,  221 ,  222 ,  223 ,  224 ,
                 -225 , -226 ,  227 ,  228 , -229 , -230 ,  231 ,  232 ,
                 -233 , -234 ,  235 ,  236 , -237 , -238 ,  239 ,  240 ,
                 -241 , -242 ,  243 ,  244 , -245 , -246 ,  247 ,  248 ,
                 -249 , -250 ,  251 ,  252 , -253 , -254 ,  255 ,  256 } ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cz62.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cz62.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cz62.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

}


/*
 * float
 */
TEST( qclab_qgates_CZ , float ) {
  test_qclab_qgates_CZ< float >() ;
}

/*
 * double
 */
TEST( qclab_qgates_CZ , double ) {
  test_qclab_qgates_CZ< double >() ;
}


/*
 * complex float
 */
TEST( qclab_qgates_CZ , complex_float ) {
  test_qclab_qgates_CZ< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_CZ , complex_double ) {
  test_qclab_qgates_CZ< std::complex< double > >() ;
}

