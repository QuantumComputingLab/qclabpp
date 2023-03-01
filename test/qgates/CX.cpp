#include <gtest/gtest.h>
#include "qclab/qgates/CX.hpp"
#include <numeric>

template <typename T>
void test_qclab_qgates_CX() {

  //
  // CXs |1> controlled
  //
  {
    qclab::qgates::CX< T >  cx ;

    EXPECT_EQ( cx.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cx.fixed() ) ;           // fixed
    EXPECT_TRUE( cx.controlled() ) ;      // controlled
    EXPECT_EQ( cx.control() , 0 ) ;       // control
    EXPECT_EQ( cx.target() , 1 ) ;        // target
    EXPECT_EQ( cx.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CX_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 0 , 1 ,
                                               0 , 0 , 1 , 0 ) ;
    EXPECT_TRUE( cx.matrix() == CX_check ) ;

    // qubit
    EXPECT_EQ( cx.qubit() , 0 ) ;

    // qubits
    auto qubits = cx.qubits() ;
    EXPECT_EQ( qubits.size() , 2 ) ;
    EXPECT_EQ( qubits[0] , 0 ) ;
    EXPECT_EQ( qubits[1] , 1 ) ;
    int qnew[] = { 5 , 3 } ;
    cx.setQubits( &qnew[0] ) ;
    EXPECT_EQ( cx.qubits()[0] , 3 ) ;
    EXPECT_EQ( cx.qubits()[1] , 5 ) ;
    qnew[0] = 0 ;
    qnew[1] = 1 ;
    cx.setQubits( &qnew[0] ) ;

    // print
    cx.print() ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cx.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cx q[0], q[1];\n" ) ;
    std::cout << qasm.str() ;

    // gate
    qclab::qgates::PauliX< T >  X ;
    EXPECT_TRUE( *cx.gate() == X ) ;
    EXPECT_TRUE( cx.gate()->matrix() == X.matrix() ) ;

    // operators == and !=
    EXPECT_TRUE(  cx != X ) ;
    EXPECT_FALSE( cx == X ) ;
    qclab::qgates::CX< T >  cx2 ;
    EXPECT_TRUE(  cx == cx2 ) ;
    EXPECT_FALSE( cx != cx2 ) ;

    // setControl, setTarget, setControlState
    cx.setControl( 3 ) ;
    EXPECT_EQ( cx.control() , 3 ) ;
    cx.setTarget( 5 ) ;
    EXPECT_EQ( cx.target() , 5 ) ;
    EXPECT_TRUE(  cx == cx2 ) ;
    EXPECT_FALSE( cx != cx2 ) ;

    cx.setControl( 4 ) ;
    EXPECT_EQ( cx.control() , 4 ) ;
    cx.setTarget( 1 ) ;
    EXPECT_EQ( cx.target() , 1 ) ;
    EXPECT_TRUE(  cx != cx2 ) ;
    EXPECT_FALSE( cx == cx2 ) ;

    qubits[0] = 1 ;
    qubits[1] = 2 ;
    cx.setQubits( &qubits[0] ) ;
    EXPECT_EQ( cx.control() , 1 ) ;
    EXPECT_EQ( cx.target() , 2 ) ;
    EXPECT_TRUE(  cx == cx2 ) ;
    EXPECT_FALSE( cx != cx2 ) ;

    cx.setControl( 0 ) ;
    EXPECT_EQ( cx.control() , 0 ) ;
    cx.setTarget( 1 ) ;
    EXPECT_EQ( cx.target() , 1 ) ;
    cx.setControlState( 0 ) ;
    EXPECT_EQ( cx.controlState() , 0 ) ;
    EXPECT_TRUE(  cx != cx2 ) ;
    EXPECT_FALSE( cx == cx2 ) ;
  }

  {
    qclab::qgates::CX< T >  cx01( 0 , 1 ) ;

    EXPECT_EQ( cx01.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cx01.fixed() ) ;           // fixed
    EXPECT_TRUE( cx01.controlled() ) ;      // controlled
    EXPECT_EQ( cx01.control() , 0 ) ;       // control
    EXPECT_EQ( cx01.target() , 1 ) ;        // target
    EXPECT_EQ( cx01.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CX_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 0 , 1 ,
                                               0 , 0 , 1 , 0 ) ;
    EXPECT_TRUE( cx01.matrix() == CX_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cx01.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cx q[0], q[1];\n" ) ;
  }

  {
    qclab::qgates::CX< T >  cx35( 3 , 5 ) ;

    EXPECT_EQ( cx35.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cx35.fixed() ) ;           // fixed
    EXPECT_TRUE( cx35.controlled() ) ;      // controlled
    EXPECT_EQ( cx35.control() , 3 ) ;       // control
    EXPECT_EQ( cx35.target() , 5 ) ;        // target
    EXPECT_EQ( cx35.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CX_check( 1 , 0 , 0 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               0 , 0 , 0 , 1 ,
                                               0 , 0 , 1 , 0 ) ;
    EXPECT_TRUE( cx35.matrix() == CX_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cx35.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cx q[3], q[5];\n" ) ;
  }

  {
    qclab::qgates::CX< T >  cx10( 1 , 0 ) ;

    EXPECT_EQ( cx10.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cx10.fixed() ) ;           // fixed
    EXPECT_TRUE( cx10.controlled() ) ;      // controlled
    EXPECT_EQ( cx10.control() , 1 ) ;       // control
    EXPECT_EQ( cx10.target() , 0 ) ;        // target
    EXPECT_EQ( cx10.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CX_check( 1 , 0 , 0 , 0 ,
                                               0 , 0 , 0 , 1 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 1 , 0 , 0 ) ;
    EXPECT_TRUE( cx10.matrix() == CX_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cx10.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cx q[1], q[0];\n" ) ;
  }

  {
    qclab::qgates::CX< T >  cx53( 5 , 3 ) ;

    EXPECT_EQ( cx53.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cx53.fixed() ) ;           // fixed
    EXPECT_TRUE( cx53.controlled() ) ;      // controlled
    EXPECT_EQ( cx53.control() , 5 ) ;       // control
    EXPECT_EQ( cx53.target() , 3 ) ;        // target
    EXPECT_EQ( cx53.controlState() , 1 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CX_check( 1 , 0 , 0 , 0 ,
                                               0 , 0 , 0 , 1 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 1 , 0 , 0 ) ;
    EXPECT_TRUE( cx53.matrix() == CX_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cx53.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "cx q[5], q[3];\n" ) ;
  }

  //
  // CXs |0> controlled
  //
  {
    qclab::qgates::CX< T >  cx01( 0 , 1 , 0 ) ;

    EXPECT_EQ( cx01.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cx01.fixed() ) ;           // fixed
    EXPECT_TRUE( cx01.controlled() ) ;      // controlled
    EXPECT_EQ( cx01.control() , 0 ) ;       // control
    EXPECT_EQ( cx01.target() , 1 ) ;        // target
    EXPECT_EQ( cx01.controlState() , 0 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CX_check( 0 , 1 , 0 , 0 ,
                                               1 , 0 , 0 , 0 ,
                                               0 , 0 , 1 , 0 ,
                                               0 , 0 , 0 , 1 ) ;
    EXPECT_TRUE( cx01.matrix() == CX_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cx01.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[0];\ncx q[0], q[1];\nx q[0];\n" ) ;
  }

  {
    qclab::qgates::CX< T >  cx10( 1 , 0 , 0 ) ;

    EXPECT_EQ( cx10.nbQubits() , 2 ) ;      // nbQubits
    EXPECT_TRUE( cx10.fixed() ) ;           // fixed
    EXPECT_TRUE( cx10.controlled() ) ;      // controlled
    EXPECT_EQ( cx10.control() , 1 ) ;       // control
    EXPECT_EQ( cx10.target() , 0 ) ;        // target
    EXPECT_EQ( cx10.controlState() , 0 ) ;  // controlState

    // matrix
    qclab::dense::SquareMatrix< T >  CX_check( 0 , 0 , 1 , 0 ,
                                               0 , 1 , 0 , 0 ,
                                               1 , 0 , 0 , 0 ,
                                               0 , 0 , 0 , 1 ) ;
    EXPECT_TRUE( cx10.matrix() == CX_check ) ;

    // toQASM
    std::stringstream qasm ;
    EXPECT_EQ( cx10.toQASM( qasm ) , 0 ) ;
    EXPECT_EQ( qasm.str() , "x q[1];\ncx q[1], q[0];\nx q[1];\n" ) ;
  }

  using V = std::vector< T > ;

  const V v2 = { 0 , 1 , 2 , 3 } ;
  const V v3 = { 0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 } ;
  V vt( 16 ) ; std::iota( vt.begin() , vt.end() , 0 ) ;
  const V v4 = vt ;
  vt = V( 32 ) ; std::iota( vt.begin() , vt.end() , 0 ) ;
  const V v5 = vt ;
  vt = V( 256 ) ; std::iota( vt.begin() , vt.end() , 0 ) ;
  const V v8 = vt ;

  //
  // CNOTs |1> controlled
  //

  // control < target
  {
    qclab::qgates::CX< T >  cx01( 0 , 1 ) ;
    qclab::qgates::CX< T >  cx12( 1 , 2 ) ;
    qclab::qgates::CX< T >  cx02( 0 , 2 ) ;
    qclab::qgates::CX< T >  cx13( 1 , 3 ) ;
    qclab::qgates::CX< T >  cx26( 2 , 6 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cx01.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 0 , 1 , 3 , 2 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cx01.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cx01.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cx01.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    cx01.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 0 , 1 , 2 , 3 , 6 , 7 , 4 , 5 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx01.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx01.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx01.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx12.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 0 , 1 , 3 , 2 , 4 , 5 , 7 , 6 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx12.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx12.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx12.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx02.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 0 , 1 , 2 , 3 , 5 , 4 , 7 , 6 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx02.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx02.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx02.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    cx12.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = {  0 ,  1 ,  2 ,  3 ,  6 ,  7 ,  4 ,  5 ,
                  8 ,  9 , 10 , 11 , 14 , 15 , 12 , 13 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cx12.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cx12.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cx12.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    cx13.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = {  0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,
                 10 , 11 ,  8 ,  9 , 14 , 15 , 12 , 13 ,
                 16 , 17 , 18 , 19 , 20 , 21 , 22 , 23 ,
                 26 , 27 , 24 , 25 , 30 , 31 , 28 , 29 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cx13.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cx13.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cx13.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    cx26.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    V check8 = {   0 ,   1 ,   2 ,   3 ,   4 ,   5 ,   6 ,   7 ,
                   8 ,   9 ,  10 ,  11 ,  12 ,  13 ,  14 ,  15 ,
                  16 ,  17 ,  18 ,  19 ,  20 ,  21 ,  22 ,  23 ,
                  24 ,  25 ,  26 ,  27 ,  28 ,  29 ,  30 ,  31 ,
                  34 ,  35 ,  32 ,  33 ,  38 ,  39 ,  36 ,  37 ,
                  42 ,  43 ,  40 ,  41 ,  46 ,  47 ,  44 ,  45 ,
                  50 ,  51 ,  48 ,  49 ,  54 ,  55 ,  52 ,  53 ,
                  58 ,  59 ,  56 ,  57 ,  62 ,  63 ,  60 ,  61 ,
                  64 ,  65 ,  66 ,  67 ,  68 ,  69 ,  70 ,  71 ,
                  72 ,  73 ,  74 ,  75 ,  76 ,  77 ,  78 ,  79 ,
                  80 ,  81 ,  82 ,  83 ,  84 ,  85 ,  86 ,  87 ,
                  88 ,  89 ,  90 ,  91 ,  92 ,  93 ,  94 ,  95 ,
                  98 ,  99 ,  96 ,  97 , 102 , 103 , 100 , 101 ,
                 106 , 107 , 104 , 105 , 110 , 111 , 108 , 109 ,
                 114 , 115 , 112 , 113 , 118 , 119 , 116 , 117 ,
                 122 , 123 , 120 , 121 , 126 , 127 , 124 , 125 ,
                 128 , 129 , 130 , 131 , 132 , 133 , 134 , 135 ,
                 136 , 137 , 138 , 139 , 140 , 141 , 142 , 143 ,
                 144 , 145 , 146 , 147 , 148 , 149 , 150 , 151 ,
                 152 , 153 , 154 , 155 , 156 , 157 , 158 , 159 ,
                 162 , 163 , 160 , 161 , 166 , 167 , 164 , 165 ,
                 170 , 171 , 168 , 169 , 174 , 175 , 172 , 173 ,
                 178 , 179 , 176 , 177 , 182 , 183 , 180 , 181 ,
                 186 , 187 , 184 , 185 , 190 , 191 , 188 , 189 ,
                 192 , 193 , 194 , 195 , 196 , 197 , 198 , 199 ,
                 200 , 201 , 202 , 203 , 204 , 205 , 206 , 207 ,
                 208 , 209 , 210 , 211 , 212 , 213 , 214 , 215 ,
                 216 , 217 , 218 , 219 , 220 , 221 , 222 , 223 ,
                 226 , 227 , 224 , 225 , 230 , 231 , 228 , 229 ,
                 234 , 235 , 232 , 233 , 238 , 239 , 236 , 237 ,
                 242 , 243 , 240 , 241 , 246 , 247 , 244 , 245 ,
                 250 , 251 , 248 , 249 , 254 , 255 , 252 , 253 } ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cx26.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cx26.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cx26.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

  // control > target
  {
    qclab::qgates::CX< T >  cx10( 1 , 0 ) ;
    qclab::qgates::CX< T >  cx21( 2 , 1 ) ;
    qclab::qgates::CX< T >  cx20( 2 , 0 ) ;
    qclab::qgates::CX< T >  cx31( 3 , 1 ) ;
    qclab::qgates::CX< T >  cx62( 6 , 2 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cx10.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 0 , 3 , 2 , 1 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cx10.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cx10.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cx10.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    cx10.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 0 , 1 , 6 , 7 , 4 , 5 , 2 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx10.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx10.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx10.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx21.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 0 , 3 , 2 , 1 , 4 , 7 , 6 , 5 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx21.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx21.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx21.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx20.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 0 , 5 , 2 , 7 , 4 , 1 , 6 , 3 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx20.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx20.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx20.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    cx21.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = {  0 ,  1 ,  6 ,  7 ,  4 ,  5 ,  2 ,  3 ,
                  8 ,  9 , 14 , 15 , 12 , 13 , 10 , 11 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cx21.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cx21.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cx21.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    cx31.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = {  0 ,  1 , 10 , 11 ,  4 ,  5 , 14 , 15 ,
                  8 ,  9 ,  2 ,  3 , 12 , 13 ,  6 ,  7 ,
                 16 , 17 , 26 , 27 , 20 , 21 , 30 , 31 ,
                 24 , 25 , 18 , 19 , 28 , 29 , 22 , 23 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cx31.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cx31.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cx31.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    cx62.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    V check8 = {   0 ,   1 ,  34 ,  35 ,   4 ,   5 ,  38 ,  39 ,
                   8 ,   9 ,  42 ,  43 ,  12 ,  13 ,  46 ,  47 ,
                  16 ,  17 ,  50 ,  51 ,  20 ,  21 ,  54 ,  55 ,
                  24 ,  25 ,  58 ,  59 ,  28 ,  29 ,  62 ,  63 ,
                  32 ,  33 ,   2 ,   3 ,  36 ,  37 ,   6 ,   7 ,
                  40 ,  41 ,  10 ,  11 ,  44 ,  45 ,  14 ,  15 ,
                  48 ,  49 ,  18 ,  19 ,  52 ,  53 ,  22 ,  23 ,
                  56 ,  57 ,  26 ,  27 ,  60 ,  61 ,  30 ,  31 ,
                  64 ,  65 ,  98 ,  99 ,  68 ,  69 , 102 , 103 ,
                  72 ,  73 , 106 , 107 ,  76 ,  77 , 110 , 111 ,
                  80 ,  81 , 114 , 115 ,  84 ,  85 , 118 , 119 ,
                  88 ,  89 , 122 , 123 ,  92 ,  93 , 126 , 127 ,
                  96 ,  97 ,  66 ,  67 , 100 , 101 ,  70 ,  71 ,
                 104 , 105 ,  74 ,  75 , 108 , 109 ,  78 ,  79 ,
                 112 , 113 ,  82 ,  83 , 116 , 117 ,  86 ,  87 ,
                 120 , 121 ,  90 ,  91 , 124 , 125 ,  94 ,  95 ,
                 128 , 129 , 162 , 163 , 132 , 133 , 166 , 167 ,
                 136 , 137 , 170 , 171 , 140 , 141 , 174 , 175 ,
                 144 , 145 , 178 , 179 , 148 , 149 , 182 , 183 ,
                 152 , 153 , 186 , 187 , 156 , 157 , 190 , 191 ,
                 160 , 161 , 130 , 131 , 164 , 165 , 134 , 135 ,
                 168 , 169 , 138 , 139 , 172 , 173 , 142 , 143 ,
                 176 , 177 , 146 , 147 , 180 , 181 , 150 , 151 ,
                 184 , 185 , 154 , 155 , 188 , 189 , 158 , 159 ,
                 192 , 193 , 226 , 227 , 196 , 197 , 230 , 231 ,
                 200 , 201 , 234 , 235 , 204 , 205 , 238 , 239 ,
                 208 , 209 , 242 , 243 , 212 , 213 , 246 , 247 ,
                 216 , 217 , 250 , 251 , 220 , 221 , 254 , 255 ,
                 224 , 225 , 194 , 195 , 228 , 229 , 198 , 199 ,
                 232 , 233 , 202 , 203 , 236 , 237 , 206 , 207 ,
                 240 , 241 , 210 , 211 , 244 , 245 , 214 , 215 ,
                 248 , 249 , 218 , 219 , 252 , 253 , 222 , 223 } ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cx62.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cx62.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cx62.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

  //
  // CNOTs |0> controlled
  //

  // control < target
  {
    qclab::qgates::CX< T >  cx01( 0 , 1 , 0 ) ;
    qclab::qgates::CX< T >  cx12( 1 , 2 , 0 ) ;
    qclab::qgates::CX< T >  cx02( 0 , 2 , 0 ) ;
    qclab::qgates::CX< T >  cx13( 1 , 3 , 0 ) ;
    qclab::qgates::CX< T >  cx26( 2 , 6 , 0 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cx01.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 1 , 0 , 2 , 3 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cx01.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cx01.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cx01.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    cx01.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 2 , 3 , 0 , 1 , 4 , 5 , 6 , 7 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx01.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx01.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx01.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx12.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 1 , 0 , 2 , 3 , 5 , 4 , 6 , 7 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx12.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx12.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx12.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx02.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 1 , 0 , 3 , 2 , 4 , 5 , 6 , 7 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx02.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx02.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx02.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    cx12.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = {  2 ,  3 ,  0 ,  1 ,  4 ,  5 ,  6 ,  7 ,
                 10 , 11 ,  8 ,  9 , 12 , 13 , 14 , 15 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cx12.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cx12.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cx12.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    cx13.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = {  2 ,  3 ,  0 ,  1 ,  6 ,  7 ,  4 ,  5 ,
                  8 ,  9 , 10 , 11 , 12 , 13 , 14 , 15 ,
                 18 , 19 , 16 , 17 , 22 , 23 , 20 , 21 ,
                 24 , 25 , 26 , 27 , 28 , 29 , 30 , 31 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cx13.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cx13.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cx13.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    cx26.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    V check8 = {   2 ,   3 ,   0 ,   1 ,   6 ,   7 ,   4 ,   5 ,
                  10 ,  11 ,   8 ,   9 ,  14 ,  15 ,  12 ,  13 ,
                  18 ,  19 ,  16 ,  17 ,  22 ,  23 ,  20 ,  21 ,
                  26 ,  27 ,  24 ,  25 ,  30 ,  31 ,  28 ,  29 ,
                  32 ,  33 ,  34 ,  35 ,  36 ,  37 ,  38 ,  39 ,
                  40 ,  41 ,  42 ,  43 ,  44 ,  45 ,  46 ,  47 ,
                  48 ,  49 ,  50 ,  51 ,  52 ,  53 ,  54 ,  55 ,
                  56 ,  57 ,  58 ,  59 ,  60 ,  61 ,  62 ,  63 ,
                  66 ,  67 ,  64 ,  65 ,  70 ,  71 ,  68 ,  69 ,
                  74 ,  75 ,  72 ,  73 ,  78 ,  79 ,  76 ,  77 ,
                  82 ,  83 ,  80 ,  81 ,  86 ,  87 ,  84 ,  85 ,
                  90 ,  91 ,  88 ,  89 ,  94 ,  95 ,  92 ,  93 ,
                  96 ,  97 ,  98 ,  99 , 100 , 101 , 102 , 103 ,
                 104 , 105 , 106 , 107 , 108 , 109 , 110 , 111 ,
                 112 , 113 , 114 , 115 , 116 , 117 , 118 , 119 ,
                 120 , 121 , 122 , 123 , 124 , 125 , 126 , 127 ,
                 130 , 131 , 128 , 129 , 134 , 135 , 132 , 133 ,
                 138 , 139 , 136 , 137 , 142 , 143 , 140 , 141 ,
                 146 , 147 , 144 , 145 , 150 , 151 , 148 , 149 ,
                 154 , 155 , 152 , 153 , 158 , 159 , 156 , 157 ,
                 160 , 161 , 162 , 163 , 164 , 165 , 166 , 167 ,
                 168 , 169 , 170 , 171 , 172 , 173 , 174 , 175 ,
                 176 , 177 , 178 , 179 , 180 , 181 , 182 , 183 ,
                 184 , 185 , 186 , 187 , 188 , 189 , 190 , 191 ,
                 194 , 195 , 192 , 193 , 198 , 199 , 196 , 197 ,
                 202 , 203 , 200 , 201 , 206 , 207 , 204 , 205 ,
                 210 , 211 , 208 , 209 , 214 , 215 , 212 , 213 ,
                 218 , 219 , 216 , 217 , 222 , 223 , 220 , 221 ,
                 224 , 225 , 226 , 227 , 228 , 229 , 230 , 231 ,
                 232 , 233 , 234 , 235 , 236 , 237 , 238 , 239 ,
                 240 , 241 , 242 , 243 , 244 , 245 , 246 , 247 ,
                 248 , 249 , 250 , 251 , 252 , 253 , 254 , 255 } ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cx26.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cx26.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cx26.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

  // control > target
  {
    qclab::qgates::CX< T >  cx10( 1 , 0 , 0 ) ;
    qclab::qgates::CX< T >  cx21( 2 , 1 , 0 ) ;
    qclab::qgates::CX< T >  cx20( 2 , 0 , 0 ) ;
    qclab::qgates::CX< T >  cx31( 3 , 1 , 0 ) ;
    qclab::qgates::CX< T >  cx62( 6 , 2 , 0 ) ;

    // nbQubits = 2
    auto vec2 = v2 ;
    cx10.apply( qclab::Op::NoTrans , 2 , vec2 ) ;
    V check2 = { 2 , 1 , 0 , 3 } ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ; T* vec2_ = vec2.data() ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cx10.apply_device( qclab::Op::NoTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    vec2 = v2 ;
    cx10.apply( qclab::Op::ConjTrans , 2 , vec2 ) ;
    EXPECT_TRUE( vec2 == check2 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec2 = v2 ;
    #pragma omp target data map(tofrom:vec2_[0:4])
    { cx10.apply_device( qclab::Op::ConjTrans , 2 , vec2_ ) ; }
    EXPECT_TRUE( vec2 == check2 ) ;
  #endif

    // nbQubits = 3
    auto vec3 = v3 ;
    cx10.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    V check3 = { 4 , 5 , 2 , 3 , 0 , 1 , 6 , 7 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ; T* vec3_ = vec3.data() ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx10.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx10.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx10.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx21.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 2 , 1 , 0 , 3 , 6 , 5 , 4 , 7 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx21.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx21.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx21.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx20.apply( qclab::Op::NoTrans , 3 , vec3 ) ;
    check3 = { 4 , 1 , 6 , 3 , 0 , 5 , 2 , 7 } ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx20.apply_device( qclab::Op::NoTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    vec3 = v3 ;
    cx20.apply( qclab::Op::ConjTrans , 3 , vec3 ) ;
    EXPECT_TRUE( vec3 == check3 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec3 = v3 ;
    #pragma omp target data map(tofrom:vec3_[0:8])
    { cx20.apply_device( qclab::Op::ConjTrans , 3 , vec3_ ) ; }
    EXPECT_TRUE( vec3 == check3 ) ;
  #endif

    // nbQubits = 4
    auto vec4 = v4 ;
    cx21.apply( qclab::Op::NoTrans , 4 , vec4 ) ;
    V check4 = {  4 ,  5 ,  2 ,  3 ,  0 ,  1 ,  6 ,  7 ,
                 12 , 13 , 10 , 11 ,  8 ,  9 , 14 , 15 } ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ; T* vec4_ = vec4.data() ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cx21.apply_device( qclab::Op::NoTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    vec4 = v4 ;
    cx21.apply( qclab::Op::ConjTrans , 4 , vec4 ) ;
    EXPECT_TRUE( vec4 == check4 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec4 = v4 ;
    #pragma omp target data map(tofrom:vec4_[0:16])
    { cx21.apply_device( qclab::Op::ConjTrans , 4 , vec4_ ) ; }
    EXPECT_TRUE( vec4 == check4 ) ;
  #endif

    // nbQubits = 5
    auto vec5 = v5 ;
    cx31.apply( qclab::Op::NoTrans , 5 , vec5 ) ;
    V check5 = {  8 ,  9 ,  2 ,  3 , 12 , 13 ,  6 ,  7 ,
                  0 ,  1 , 10 , 11 ,  4 ,  5 , 14 , 15 ,
                 24 , 25 , 18 , 19 , 28 , 29 , 22 , 23 ,
                 16 , 17 , 26 , 27 , 20 , 21 , 30 , 31 } ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ; T* vec5_ = vec5.data() ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cx31.apply_device( qclab::Op::NoTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    vec5 = v5 ;
    cx31.apply( qclab::Op::ConjTrans , 5 , vec5 ) ;
    EXPECT_TRUE( vec5 == check5 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec5 = v5 ;
    #pragma omp target data map(tofrom:vec5_[0:32])
    { cx31.apply_device( qclab::Op::ConjTrans , 5 , vec5_ ) ; }
    EXPECT_TRUE( vec5 == check5 ) ;
  #endif

    // nbQubits = 8
    auto vec8 = v8 ;
    cx62.apply( qclab::Op::NoTrans , 8 , vec8 ) ;
    V check8 = {  32 ,  33 ,   2 ,   3 ,  36 ,  37 ,   6 ,   7 ,
                  40 ,  41 ,  10 ,  11 ,  44 ,  45 ,  14 ,  15 ,
                  48 ,  49 ,  18 ,  19 ,  52 ,  53 ,  22 ,  23 ,
                  56 ,  57 ,  26 ,  27 ,  60 ,  61 ,  30 ,  31 ,
                   0 ,   1 ,  34 ,  35 ,   4 ,   5 ,  38 ,  39 ,
                   8 ,   9 ,  42 ,  43 ,  12 ,  13 ,  46 ,  47 ,
                  16 ,  17 ,  50 ,  51 ,  20 ,  21 ,  54 ,  55 ,
                  24 ,  25 ,  58 ,  59 ,  28 ,  29 ,  62 ,  63 ,
                  96 ,  97 ,  66 ,  67 , 100 , 101 ,  70 ,  71 ,
                 104 , 105 ,  74 ,  75 , 108 , 109 ,  78 ,  79 ,
                 112 , 113 ,  82 ,  83 , 116 , 117 ,  86 ,  87 ,
                 120 , 121 ,  90 ,  91 , 124 , 125 ,  94 ,  95 ,
                  64 ,  65 ,  98 ,  99 ,  68 ,  69 , 102 , 103 ,
                  72 ,  73 , 106 , 107 ,  76 ,  77 , 110 , 111 ,
                  80 ,  81 , 114 , 115 ,  84 ,  85 , 118 , 119 ,
                  88 ,  89 , 122 , 123 ,  92 ,  93 , 126 , 127 ,
                 160 , 161 , 130 , 131 , 164 , 165 , 134 , 135 ,
                 168 , 169 , 138 , 139 , 172 , 173 , 142 , 143 ,
                 176 , 177 , 146 , 147 , 180 , 181 , 150 , 151 ,
                 184 , 185 , 154 , 155 , 188 , 189 , 158 , 159 ,
                 128 , 129 , 162 , 163 , 132 , 133 , 166 , 167 ,
                 136 , 137 , 170 , 171 , 140 , 141 , 174 , 175 ,
                 144 , 145 , 178 , 179 , 148 , 149 , 182 , 183 ,
                 152 , 153 , 186 , 187 , 156 , 157 , 190 , 191 ,
                 224 , 225 , 194 , 195 , 228 , 229 , 198 , 199 ,
                 232 , 233 , 202 , 203 , 236 , 237 , 206 , 207 ,
                 240 , 241 , 210 , 211 , 244 , 245 , 214 , 215 ,
                 248 , 249 , 218 , 219 , 252 , 253 , 222 , 223 ,
                 192 , 193 , 226 , 227 , 196 , 197 , 230 , 231 ,
                 200 , 201 , 234 , 235 , 204 , 205 , 238 , 239 ,
                 208 , 209 , 242 , 243 , 212 , 213 , 246 , 247 ,
                 216 , 217 , 250 , 251 , 220 , 221 , 254 , 255 } ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ; T* vec8_ = vec8.data() ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cx62.apply_device( qclab::Op::NoTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif

    vec8 = v8 ;
    cx62.apply( qclab::Op::ConjTrans , 8 , vec8 ) ;
    EXPECT_TRUE( vec8 == check8 ) ;
  #ifdef QCLAB_OMP_OFFLOADING
    vec8 = v8 ;
    #pragma omp target data map(tofrom:vec8_[0:256])
    { cx62.apply_device( qclab::Op::ConjTrans , 8 , vec8_ ) ; }
    EXPECT_TRUE( vec8 == check8 ) ;
  #endif
  }

}


/*
 * float
 */
TEST( qclab_qgates_CX , float ) {
  test_qclab_qgates_CX< float >() ;
}

/*
 * double
 */
TEST( qclab_qgates_CX , double ) {
  test_qclab_qgates_CX< double >() ;
}


/*
 * complex float
 */
TEST( qclab_qgates_CX , complex_float ) {
  test_qclab_qgates_CX< std::complex< float > >() ;
}

/*
 * complex double
 */
TEST( qclab_qgates_CX , complex_double ) {
  test_qclab_qgates_CX< std::complex< double > >() ;
}

