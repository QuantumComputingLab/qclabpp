#include "qclab/io/QASMFile.hpp"
#include <algorithm>
#include <cassert>

namespace qclab::io {

  // QASMFile
  QASMFile::QASMFile( const std::string filename )
  : stream_( std::ifstream( filename ) )
  {
    assert( stream_.good() ) ;
    parse() ;
  } // QASMFile(filename)

  // parse
  void QASMFile::parse() {

    bool hasQreg = false ;

    // loop over lines in file
    while ( not stream_.eof() ) {

      // new line + remove all whitespaces
      std::string str ;
      std::getline( stream_ , str ) ;
      str.erase( remove_if( str.begin() , str.end() , isspace ) , str.end() ) ;

      // skip blank and comment lines
      if ( str.length() < 1 ) { continue ; }
      if ( str.substr( 0 , 2 ) == "//" ) { continue ; }

      // split commands on this line
      auto commands = split< std::string >( str , ";" ) ;

      // parse commands
      for ( auto& command : commands ) {
        if ( !hasQreg ) {
          if ( command.substr( 0 , 4 ) == "qreg" ) {
            // set quantum register
            command = command.substr( 4 ) ;
            qregName_ = left_of( command , "[" ) ;
            nbQubits_ = read_value< int >( left_of( command , "]" ) ) ;
            hasQreg = true ;
          }
        } else {
          // read gates
          gateCommands_.push_back( command ) ;
        }
      }

    }

  } // parse()

} // namespace qclab::io

