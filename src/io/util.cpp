#include "qclab/io/util.hpp"

namespace qclab::io {

  /// Reads int from `str`.
  template <>
  int         read_value( const std::string& str ) { return std::stoi( str ) ; }

  /// Reads float from `str`.
  template <>
  float       read_value( const std::string& str ) { return std::stof( str ) ; }

  /// Reads double from `str`.
  template <>
  double      read_value( const std::string& str ) { return std::stod( str ) ; }

  /// Reads string from `str`.
  template <>
  std::string read_value( const std::string& str ) { return str ; }

} // namespace qclab::io

