//  (C) Copyright Roel Van Beeumen and Daan Camps 2022.

#pragma once

#include <string>
#include <vector>

namespace qclab::io {

  /// Reads value for `str`.
  template <typename T>
  T read_value( const std::string& str ) ;

  /// Parses the string `str` left of the delimiter `delimiter`.
  inline std::string left_of( std::string& str , const std::string delimiter ) {
    auto pos = str.find_first_of( delimiter ) ;
    if ( pos == 0 ) {
      return "" ;
    } else if ( pos == std::string::npos ) {
      return str ;
    }
    std::string left_of = str.substr( 0 , pos ) ;
    str = str.substr( pos + 1 ) ;
    return left_of ;
  }

  /// Splits the string `str` in a vector using the delimiter `delimiter`.
  template <typename T>
  std::vector< T > split( std::string& str , const std::string delimiter ) {
    std::vector< T > vec ;
    while ( str.length() > 0 ) {
      vec.push_back( read_value< T >( left_of( str , delimiter ) ) ) ;
    }
    return vec ;
  }

} // namespace qclab::io

