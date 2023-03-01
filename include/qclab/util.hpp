//  (C) Copyright Roel Van Beeumen 2021.

#pragma once

#include <cassert>
#include <complex>
#include <iostream>

namespace qclab {

  /**
   * Operation to perform.
   * \ingroup Enumerations
   */
  enum class Op : char
  {
    NoTrans   = 'N' ,  ///< No transpose
    Trans     = 'T' ,  ///< Transpose
    ConjTrans = 'C'    ///< Conjugate transpose
  } ; // enum class Op

  /**
   * Which side to apply the operation on.
   * \ingroup Enumerations
   */
  enum class Side : char
  {
    Left  = 'L' ,  ///< Left
    Right = 'R'    ///< Right
  } ; // enum class Side


  /// Default helper class for real_t.
  template <typename T>
  struct real_t_helper {
    typedef T real_type ;  ///< Real type
  } ;

  /// Template specialized helper class for real_t.
  template <typename T>
  struct real_t_helper< std::complex< T > > {
    typedef T real_type ;  ///< Real type
  } ;

  /// Typedef for real type.
  template <typename T>
  using real_t = typename real_t_helper< T >::real_type ;


  /// Default helper class for is_complex.
  template <typename T>
  struct is_complex_helper : std::false_type { } ;

  /// Template specialized helper class for is_complex.
  template <typename T>
  struct is_complex_helper< std::complex< T > > : std::true_type { } ;

  /// Returns type of is_complex.
  template <typename T>
  struct is_complex : is_complex_helper< T >::type { } ;

  /// Checks if T is complex.
  template <typename T>
  inline constexpr bool is_complex_v = is_complex< T >::value ;


  /// Print real value.
  template <typename T>
  void print( const T val ) {
    if ( val == T(0) ) {
      std::printf( "         0" ) ;
    } else {
      std::printf( "%10.4f" , val ) ;
    }
  } // print

  /// Print complex value.
  template <typename T>
  void print( const std::complex< T > val ) {
    const T real = std::real( val ) ;
    const T imag = std::abs( std::imag( val ) ) ;
    const char sign = std::imag( val ) >= T(0) ? '+' : '-' ;
    if ( std::abs( val ) == T(0) ) {
      std::printf( "                   0" ) ;
    } else if ( std::abs( real ) == T(0) ) {
      std::printf( "         0 %c %6.4fi" , sign , imag ) ;
    } else if ( imag == T(0) ) {
      std::printf( "%10.4f %c      0i" , real , sign ) ;
    } else {
      std::printf( "%10.4f %c %6.4fi" , real , sign , imag ) ;
    }
  } // print

  /// Print 2x2 matrix.
  template <typename T>
  void printMatrix2x2( const T& mat ) {
    assert( mat.rows() == 2 ) ;
    assert( mat.cols() == 2 ) ;
    print( mat(0,0) ) ; print( mat(0,1) ) ; std::cout << std::endl ;
    print( mat(1,0) ) ; print( mat(1,1) ) ; std::cout << std::endl ;
  } // printMatrix2x2

  /// Print 4x4 matrix.
  template <typename T>
  void printMatrix4x4( const T& mat ) {
    assert( mat.rows() == 4 ) ;
    assert( mat.cols() == 4 ) ;
    print( mat(0,0) ); print( mat(0,1) ); print( mat(0,2) ); print( mat(0,3) );
    std::cout << std::endl ;
    print( mat(1,0) ); print( mat(1,1) ); print( mat(1,2) ); print( mat(1,3) );
    std::cout << std::endl ;
    print( mat(2,0) ); print( mat(2,1) ); print( mat(2,2) ); print( mat(2,3) );
    std::cout << std::endl ;
    print( mat(3,0) ); print( mat(3,1) ); print( mat(3,2) ); print( mat(3,3) );
    std::cout << std::endl ;
  } // printMatrix4x4

  /// Print matrix.
  template <typename T>
  void printMatrix( const T& mat ) {
   for ( size_t r = 0; r < mat.rows(); ++r ) {
     for ( size_t c = 0; c < mat.cols(); ++c ) {
       print( mat(r,c) ) ;
     }
     std::cout << std::endl ;
   }
 }

} // namespace qclab

