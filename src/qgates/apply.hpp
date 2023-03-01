//  (C) Copyright Roel Van Beeumen and Daan Camps 2022.

#pragma once

#include <tuple>

namespace qclab::qgates {

  // lambda_QGate1
  template <typename T>
  auto lambda_QGate1( Op op , qclab::dense::SquareMatrix< T > mat1 ,
                      T* vector ) {
    assert( mat1.size() == 2 ) ;
    // operation
    qclab::dense::operateInPlace( op , mat1 ) ;
    const T m11 = mat1( 0 , 0 ) ; const T m12 = mat1( 0 , 1 ) ;
    const T m21 = mat1( 1 , 0 ) ; const T m22 = mat1( 1 , 1 ) ;
    // matvec
    auto f = [=] ( const uint64_t a , const uint64_t b ) {
      const T x1 = vector[a] ;
      const T x2 = vector[b] ;
      vector[a] = m11 * x1 + m12 * x2 ;
      vector[b] = m21 * x1 + m22 * x2 ;
    } ;
    return f ;
  }

  // lambda_Hadamard
  template <typename T>
  auto lambda_Hadamard( Op op , T* vector ) {
    // operation
    using R = qclab::real_t< T > ;
    const R sqrt2 = R(1) / std::sqrt( R(2) ) ;
    // matvec
    auto f = [=] ( const uint64_t a , const uint64_t b ) {
      const T x1 = vector[a] ;
      const T x2 = vector[b] ;
      vector[a] = sqrt2 * x1 + sqrt2 * x2 ;
      vector[b] = sqrt2 * x1 - sqrt2 * x2 ;
    } ;
    return f ;
  }

  // lambda_PauliX
  template <typename T>
  auto lambda_PauliX( Op op , T* vector ) {
    // matvec
    auto f = [=] ( const uint64_t a , const uint64_t b ) {
      std::swap( vector[a] , vector[b] ) ;
    } ;
    return f ;
  }

  // lambda_PauliY
  template <typename T>
  auto lambda_PauliY( Op op , T* vector ) {
    // operation
    T lambda1 = T( 0 , -1 ) ;
    T lambda2 = T( 0 ,  1 ) ;
    if ( op == Op::Trans ) std::swap( lambda1 , lambda2 ) ;
    // matvec
    auto f = [=] ( const uint64_t a , const uint64_t b ) {
      std::swap( vector[a] , vector[b] ) ;
      vector[a] *= lambda1 ;
      vector[b] *= lambda2 ;
    } ;
    return f ;
  }

  // lambda_PauliZ
  template <typename T>
  auto lambda_PauliZ( Op op , T* vector ) {
    // matvec
    auto f = [=] ( const uint64_t b ) {
      vector[b] = -vector[b] ;
    } ;
    return f ;
  }

  // lambda_RotationX
  template <typename T, typename R = qclab::real_t< T >>
  auto lambda_RotationX( Op op , const R cos , const R sin , T* vector ) {
    // operation
    const R c = cos ;
    T s = T( 0 , -sin ) ;
    if ( op == Op::ConjTrans ) s = std::conj( s ) ;
    // matvec
    auto f = [=] ( const uint64_t a , const uint64_t b ) {
      const T x1 = vector[a] ;
      const T x2 = vector[b] ;
      vector[a] = c * x1 + s * x2 ;
      vector[b] = s * x1 + c * x2 ;
    } ;
    return f ;
  }

  // lambda_RotationY
  template <typename T, typename R = qclab::real_t< T >>
  auto lambda_RotationY( Op op , const R cos , const R sin , T* vector ) {
    // operation
    const R c = cos ;
    R s1 = -sin ;
    R s2 =  sin ;
    if ( op != Op::NoTrans ) std::swap( s1 , s2 ) ;
    // matvec
    auto f = [=] ( const uint64_t a , const uint64_t b ) {
      const T x1 = vector[a] ;
      const T x2 = vector[b] ;
      vector[a] =  c * x1 + s1 * x2 ;
      vector[b] = s2 * x1 +  c * x2 ;
    } ;
    return f ;
  }

  // lambda_RotationZ
  template <typename T, typename R = qclab::real_t< T >>
  auto lambda_RotationZ( Op op , const R cos , const R sin , T* vector ) {
    // operation
    T lambda1 = T( cos , -sin ) ;
    T lambda2 = T( cos ,  sin ) ;
    if ( op == Op::ConjTrans ) std::swap( lambda1 , lambda2 ) ;
    // matvec
    auto f = [=] ( const uint64_t a , const uint64_t b ) {
      vector[a] *= lambda1 ;
      vector[b] *= lambda2 ;
    } ;
    return f ;
  }

  // lambda_Phase
  template <typename T>
  auto lambda_Phase( Op op , T& lambda , T* vector ) {
    // operation
    if ( op == Op::ConjTrans ) lambda = std::conj( lambda ) ;
    // matvec
    auto f = [=] ( const uint64_t b ) {
      vector[b] *= lambda ;
    } ;
    return f ;
  }

  // lambda_Phase45
  template <typename T>
  auto lambda_Phase45( Op op , T* vector ) {
    // operation
    using R = qclab::real_t< T > ;
    const R sqrt2 = R(1) / std::sqrt( R(2) ) ;
    T lambda = T( sqrt2 , sqrt2 ) ;
    if ( op == Op::ConjTrans ) lambda = std::conj( lambda ) ;
    // matvec
    auto f = [=] ( const uint64_t b ) {
     vector[b] *= lambda ;
    } ;
    return f ;
  }

  // lambda_Phase90
  template <typename T>
  auto lambda_Phase90( Op op , T* vector ) {
    // operation
    T lambda = T( 0 , 1 ) ;
    if ( op == Op::ConjTrans ) lambda = std::conj( lambda ) ;
    // matvec
    auto f = [=] ( const uint64_t b ) {
      vector[b] *= lambda ;
    } ;
    return f ;
  }

  // lambda_QGate2
  template <typename T>
  auto lambda_QGate2( Op op , qclab::dense::SquareMatrix< T > mat2 ,
                      T* vector ) {
    assert( mat2.size() == 4 ) ;
    // operation
    qclab::dense::operateInPlace( op , mat2 ) ;
    const T m11 = mat2( 0 , 0 ) ; const T m12 = mat2( 0 , 1 ) ;
    const T m13 = mat2( 0 , 2 ) ; const T m14 = mat2( 0 , 3 ) ;
    const T m21 = mat2( 1 , 0 ) ; const T m22 = mat2( 1 , 1 ) ;
    const T m23 = mat2( 1 , 2 ) ; const T m24 = mat2( 1 , 3 ) ;
    const T m31 = mat2( 2 , 0 ) ; const T m32 = mat2( 2 , 1 ) ;
    const T m33 = mat2( 2 , 2 ) ; const T m34 = mat2( 2 , 3 ) ;
    const T m41 = mat2( 3 , 0 ) ; const T m42 = mat2( 3 , 1 ) ;
    const T m43 = mat2( 3 , 2 ) ; const T m44 = mat2( 3 , 3 ) ;
    // matvec
    auto f = [=] ( const uint64_t a , const uint64_t b ,
                   const uint64_t c , const uint64_t d ) {
      const T x1 = vector[a] ;
      const T x2 = vector[b] ;
      const T x3 = vector[c] ;
      const T x4 = vector[d] ;
      vector[a] = m11 * x1 + m12 * x2 + m13 * x3 + m14 * x4 ;
      vector[b] = m21 * x1 + m22 * x2 + m23 * x3 + m24 * x4 ;
      vector[c] = m31 * x1 + m32 * x2 + m33 * x3 + m34 * x4 ;
      vector[d] = m41 * x1 + m42 * x2 + m43 * x3 + m44 * x4 ;
    } ;
    return f ;
  }

  // lambda_SWAP
  template <typename T>
  auto lambda_SWAP( Op op , T* vector ) {
    // matvec
    auto f = [=] ( const uint64_t b , const uint64_t c ) {
      std::swap( vector[b] , vector[c] ) ;
    } ;
    return f ;
  }

  // lambda_iSWAP
  template <typename T>
  auto lambda_iSWAP( Op op , T* vector ) {
    // operation
    T lambda = T( 0 , 1 ) ;
    if ( op == Op::ConjTrans ) lambda = T( 0 , -1 ) ;
    // matvec
    auto f = [=] ( const uint64_t b , const uint64_t c ) {
      std::swap( vector[b] , vector[c] ) ;
      vector[b] *= lambda ;
      vector[c] *= lambda ;
    } ;
    return f ;
  }

  // masks for 1-qubit gates
  inline
  std::tuple< uint64_t , uint64_t > masks( const int nbQubits ,
                                           const int qubit ) {
    assert( nbQubits >= 1 ) ;
    assert( qubit < nbQubits ) ;
    // mask lengths
    const int nL = qubit ;
    const int nR = nbQubits - qubit - 1 ;
    // masks
    const uint64_t mL = ( ( 1ULL << nL ) - 1 ) << nR ;
    const uint64_t mR =   ( 1ULL << nR ) - 1 ;
    return { mL , mR } ;
  }

  // masks for 2-qubit gates
  inline
  std::tuple< uint64_t , uint64_t , uint64_t > masks( const int nbQubits ,
                                                      const int qubit0 ,
                                                      const int qubit1 ) {
    assert( nbQubits >= 2 ) ;
    assert( qubit0 < nbQubits ) ; assert( qubit1 < nbQubits ) ;
    // mask lengths
    const int nL = qubit0 ;
    const int nC = qubit1 - qubit0 - 1 ;
    const int nR = nbQubits - qubit1 - 1 ;
    // masks
    const uint64_t mL = ( ( 1ULL << nL ) - 1 ) << ( nR + nC ) ;
    const uint64_t mC = ( ( 1ULL << nC ) - 1 ) <<   nR ;
    const uint64_t mR =   ( 1ULL << nR ) - 1 ;
    return { mL , mC , mR } ;
  }

  template <typename F>
  void apply2( const int nbQubits , const int qubit , F& lambda ) {
    // masks
    const auto [ mL , mR ] = masks( nbQubits , qubit ) ;
    // indices
    const uint64_t n  = 1ULL << ( nbQubits - 1 ) ;
    const uint64_t b1 = 1ULL << ( nbQubits - qubit - 1 ) ;
    // matvec
    #pragma omp parallel for
    for ( uint64_t k = 0; k < n; k++ ) {
      const uint64_t a = ( (k & mL) << 1 ) | (k & mR) ;
      const uint64_t b = a | b1 ;
      lambda( a , b ) ;
    }
  }

  template <typename F>
  void apply2b( const int nbQubits , const int qubit , F& lambda ) {
    // masks
    const auto [ mL , mR ] = masks( nbQubits , qubit ) ;
    // indices
    const uint64_t n  = 1ULL << ( nbQubits - 1 ) ;
    const uint64_t b1 = 1ULL << ( nbQubits - qubit - 1 ) ;
    // matvec
    #pragma omp parallel for
    for ( uint64_t k = 0; k < n; k++ ) {
      const uint64_t b = ( (k & mL) << 1 ) | b1 | (k & mR) ;
      lambda( b ) ;
    }
  }

  template <typename F>
  void apply4( const int nbQubits , const int qubit0 , const int qubit1 ,
               F& lambda ) {
    // masks
    const auto [ mL , mC , mR ] = masks( nbQubits , qubit0 , qubit1 ) ;
    // indices
    const uint64_t  n = 1ULL << ( nbQubits - 2 ) ;
    const uint64_t b1 = 1ULL << ( nbQubits - qubit0 - 1 ) ;
    const uint64_t c1 = 1ULL << ( nbQubits - qubit1 - 1 ) ;
    const uint64_t d1 = b1 | c1 ;
    // matvec
    #pragma omp parallel for
    for ( uint64_t k = 0; k < n; k++ ) {
      const uint64_t a = (k & mR) + ( (k & mC) << 1 ) + ( (k & mL) << 2 ) ;
      const uint64_t b = a | b1 ;
      const uint64_t c = a | c1 ;
      const uint64_t d = a | d1 ;
      lambda( a , b , c , d ) ;
    }
  }

  template <typename F>
  void apply4( const int nbQubits , const int qubit0 , const int qubit1 ,
               const int control , const int target , const int controlState ,
               F& lambda ) {
    // masks
    const auto [ mL , mC , mR ] = masks( nbQubits , qubit0 , qubit1 ) ;
    // control / target
    const uint64_t n  = 1ULL << ( nbQubits - 2 ) ;
          uint64_t a1 = 0ULL ;
          uint64_t b1 = 1ULL << ( nbQubits - target - 1 ) ;
    if ( controlState == 1 ) {
      a1 += 1ULL << ( nbQubits - control - 1 ) ;
      b1 += 1ULL << ( nbQubits - control - 1 ) ;
    }
    // matvec
    #pragma omp parallel for
    for ( uint64_t k = 0; k < n; k++ ) {
      const uint64_t i = ( (k & mL) << 2 ) | ( (k & mC) << 1 ) | (k & mR) ;
      const uint64_t a = i | a1 ;
      const uint64_t b = i | b1 ;
      lambda( a , b ) ;
    }
  }

  template <typename F>
  void apply4b( const int nbQubits , const int qubit0 , const int qubit1 ,
                const int control , const int target , const int controlState ,
                F& lambda ) {
    // masks
    const auto [ mL , mC , mR ] = masks( nbQubits , qubit0 , qubit1 ) ;
    // control / target
    const uint64_t n  = 1ULL << ( nbQubits - 2 ) ;
          uint64_t b1 = 1ULL << ( nbQubits - target - 1 ) ;
    if ( controlState == 1 ) {
      b1 += 1ULL << ( nbQubits - control - 1 ) ;
    }
    // matvec
    #pragma omp parallel for
    for ( uint64_t k = 0; k < n; k++ ) {
      const uint64_t b = ( (k & mL) << 2 ) | ( (k & mC) << 1 ) | (k & mR) | b1 ;
      lambda( b ) ;
    }
  }

  template <typename F>
  void apply4bc( const int nbQubits , const int qubit0 , const int qubit1 ,
                 F& lambda ) {
    // masks
    const auto [ mL , mC , mR ] = masks( nbQubits , qubit0 , qubit1 ) ;
    // indices
    const uint64_t n  = 1ULL << ( nbQubits - 2 ) ;
    const uint64_t b1 = 1ULL << ( nbQubits - qubit0 - 1 ) ;
    const uint64_t c1 = 1ULL << ( nbQubits - qubit1 - 1 ) ;
    // matvec
    #pragma omp parallel for
    for ( uint64_t k = 0; k < n; k++ ) {
      const uint64_t a = ( (k & mL) << 2 ) | ( (k & mC) << 1 ) | (k & mR) ;
      const uint64_t b = a | b1 ;
      const uint64_t c = a | c1 ;
      lambda( b , c ) ;
    }
  }

#ifdef QCLAB_OMP_OFFLOADING
  template <typename F>
  void apply_device2( const int nbQubits , const int qubit , F& lambda ) {
    assert( nbQubits >= 1 ) ;
    assert( qubit < nbQubits ) ;
    // mask lengths
    const int nL = qubit ;
    const int nR = nbQubits - qubit - 1 ;
    // masks
    const uint64_t mL = ( ( 1ULL << nL ) - 1 ) << nR ;
    const uint64_t mR =   ( 1ULL << nR ) - 1 ;
    // indices
    const uint64_t n  = 1ULL << ( nbQubits - 1 ) ;
    const uint64_t b1 = 1ULL << ( nbQubits - qubit - 1 ) ;
    // matvec
    #pragma omp target teams distribute parallel for
    for ( uint64_t k = 0; k < n; k++ ) {
      const uint64_t a = ( (k & mL) << 1 ) | (k & mR) ;
      const uint64_t b = a | b1 ;
      lambda( a , b ) ;
    }
  }

  template <typename F>
  void apply_device2b( const int nbQubits , const int qubit , F& lambda ) {
    assert( nbQubits >= 1 ) ;
    assert( qubit < nbQubits ) ;
    // mask lengths
    const int nL = qubit ;
    const int nR = nbQubits - qubit - 1 ;
    // masks
    const uint64_t mL = ( ( 1ULL << nL ) - 1 ) << nR ;
    const uint64_t mR =   ( 1ULL << nR ) - 1 ;
    // indices
    const uint64_t n  = 1ULL << ( nbQubits - 1 ) ;
    const uint64_t b1 = 1ULL << ( nbQubits - qubit - 1 ) ;
    // matvec
    #pragma omp target teams distribute parallel for
    for ( uint64_t k = 0; k < n; k++ ) {
      const uint64_t b = ( (k & mL) << 1 ) | b1 | (k & mR) ;
      lambda( b ) ;
    }
  }

  template <typename F>
  void apply_device4( const int nbQubits , const int qubit0 ,
                      const int qubit1 , F& lambda ) {
    assert( nbQubits >= 2 ) ;
    assert( qubit0 < nbQubits ) ; assert( qubit1 < nbQubits ) ;
    // mask lengths
    const int nL = qubit0 ;
    const int nC = qubit1 - qubit0 - 1 ;
    const int nR = nbQubits - qubit1 - 1 ;
    // masks
    const uint64_t mL = ( ( 1ULL << nL ) - 1 ) << ( nbQubits - qubit0 - 2 ) ;
    const uint64_t mC = ( ( 1ULL << nC ) - 1 ) << ( nbQubits - qubit1 - 1 ) ;
    const uint64_t mR =   ( 1ULL << nR ) - 1 ;
    // indices
    const uint64_t  n = 1ULL << ( nbQubits - 2 ) ;
    const uint64_t b1 = 1ULL << ( nbQubits - qubit0 - 1 ) ;
    const uint64_t c1 = 1ULL << ( nbQubits - qubit1 - 1 ) ;
    const uint64_t d1 = b1 | c1 ;
    // matvec
    #pragma omp target teams distribute parallel for
    for ( uint64_t k = 0; k < n; k++ ) {
      const uint64_t a = (k & mR) + ( (k & mC) << 1 ) + ( (k & mL) << 2 ) ;
      const uint64_t b = a | b1 ;
      const uint64_t c = a | c1 ;
      const uint64_t d = a | d1 ;
      lambda( a , b , c , d ) ;
    }
  }

  template <typename F>
  void apply_device4( const int nbQubits , const int qubit0 ,
                      const int qubit1 , const int control ,
                      const int target , const int controlState , F& lambda ) {
    assert( nbQubits >= 2 ) ;
    assert( qubit0 < nbQubits ) ; assert( qubit1 < nbQubits ) ;
    // mask lengths
    const int nL = qubit0 ;
    const int nC = qubit1 - qubit0 - 1 ;
    const int nR = nbQubits - qubit1 - 1 ;
    // masks
    const uint64_t mL = ( ( 1ULL << nL ) - 1 ) << ( nbQubits - qubit0 - 2 ) ;
    const uint64_t mC = ( ( 1ULL << nC ) - 1 ) << ( nbQubits - qubit1 - 1 ) ;
    const uint64_t mR =   ( 1ULL << nR ) - 1 ;
    // control / target
    const uint64_t n  = 1ULL << ( nbQubits - 2 ) ;
          uint64_t a1 = 0ULL ;
          uint64_t b1 = 1ULL << ( nbQubits - target - 1 ) ;
    if ( controlState == 1 ) {
      a1 += 1ULL << ( nbQubits - control - 1 ) ;
      b1 += 1ULL << ( nbQubits - control - 1 ) ;
    }
    // matvec
    #pragma omp target teams distribute parallel for
    for ( uint64_t k = 0; k < n; k++ ) {
      const uint64_t i = ( (k & mL) << 2 ) | ( (k & mC) << 1 ) | (k & mR) ;
      const uint64_t a = i | a1 ;
      const uint64_t b = i | b1 ;
      lambda( a , b ) ;
    }
  }

  template <typename F>
  void apply_device4b( const int nbQubits , const int qubit0 ,
                       const int qubit1 , const int control ,
                       const int target , const int controlState , F& lambda ) {
    assert( nbQubits >= 2 ) ;
    assert( qubit0 < nbQubits ) ; assert( qubit1 < nbQubits ) ;
    // mask lengths
    const int nL = qubit0 ;
    const int nC = qubit1 - qubit0 - 1 ;
    const int nR = nbQubits - qubit1 - 1 ;
    // masks
    const uint64_t mL = ( ( 1ULL << nL ) - 1 ) << ( nbQubits - qubit0 - 2 ) ;
    const uint64_t mC = ( ( 1ULL << nC ) - 1 ) << ( nbQubits - qubit1 - 1 ) ;
    const uint64_t mR =   ( 1ULL << nR ) - 1 ;
    // control / target
    const uint64_t n  = 1ULL << ( nbQubits - 2 ) ;
          uint64_t b1 = 1ULL << ( nbQubits - target - 1 ) ;
    if ( controlState == 1 ) {
      b1 += 1ULL << ( nbQubits - control - 1 ) ;
    }
    // matvec
    #pragma omp target teams distribute parallel for
    for ( uint64_t k = 0; k < n; k++ ) {
      const uint64_t b = ( (k & mL) << 2 ) | ( (k & mC) << 1 ) | (k & mR) | b1 ;
      lambda( b ) ;
    }
  }

  template <typename F>
  void apply_device4bc( const int nbQubits , const int qubit0 ,
                        const int qubit1 , F& lambda ) {
    assert( nbQubits >= 2 ) ;
    assert( qubit0 < nbQubits ) ; assert( qubit1 < nbQubits ) ;
    // mask lengths
    const int nL = qubit0 ;
    const int nC = qubit1 - qubit0 - 1 ;
    const int nR = nbQubits - qubit1 - 1 ;
    // masks
    const uint64_t mL = ( ( 1ULL << nL ) - 1 ) << ( nbQubits - qubit0 - 2 ) ;
    const uint64_t mC = ( ( 1ULL << nC ) - 1 ) << ( nbQubits - qubit1 - 1 ) ;
    const uint64_t mR =   ( 1ULL << nR ) - 1 ;
    // indices
    const uint64_t n  = 1ULL << ( nbQubits - 2 ) ;
    const uint64_t b1 = 1ULL << ( nbQubits - qubit0 - 1 ) ;
    const uint64_t c1 = 1ULL << ( nbQubits - qubit1 - 1 ) ;
    // matvec
    #pragma omp target teams distribute parallel for
    for ( uint64_t k = 0; k < n; k++ ) {
      const uint64_t a = ( (k & mL) << 2 ) | ( (k & mC) << 1 ) | (k & mR) ;
      const uint64_t b = a | b1 ;
      const uint64_t c = a | c1 ;
      lambda( b , c ) ;
    }
  }
#endif

} // namespace qclab::qgates

