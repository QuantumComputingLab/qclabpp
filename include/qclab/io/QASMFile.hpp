//  (C) Copyright Roel Van Beeumen and Daan Camps 2022.

#pragma once

#include "qclab/io/util.hpp"
#include "qclab/qgates/Hadamard.hpp"
#include "qclab/qgates/PauliX.hpp"
#include "qclab/qgates/PauliY.hpp"
#include "qclab/qgates/PauliZ.hpp"
#include "qclab/qgates/RotationX.hpp"
#include "qclab/qgates/RotationY.hpp"
#include "qclab/qgates/RotationZ.hpp"
#include "qclab/qgates/Phase.hpp"
#include "qclab/qgates/CX.hpp"
#include "qclab/qgates/CY.hpp"
#include "qclab/qgates/CZ.hpp"
#include "qclab/qgates/SWAP.hpp"
#include <fstream>
#include <memory>
#include <iostream>

/**
 * Namespace qclab::io.
 */
namespace qclab::io {

  /**
   * \class QASMFile
   * \brief Class for parsing qasm files.
   */
  class QASMFile
  {

    public:
      /// Constructs a QASM file parser for the given filename.
      QASMFile( const std::string filename ) ;

      /// Returns the number of qubits.
      int nbQubits() const { return nbQubits_ ; }

      /// Returns the gates.
      template <typename G>
      std::vector< std::unique_ptr< G > > gates() const {

        // types
        using T = typename G::value_type ;
        using H = qclab::qgates::Hadamard< T > ;
        using X = qclab::qgates::PauliX< T > ;
        using Y = qclab::qgates::PauliY< T > ;
        using Z = qclab::qgates::PauliZ< T > ;
        using RX = qclab::qgates::RotationX< T > ;
        using RY = qclab::qgates::RotationY< T > ;
        using RZ = qclab::qgates::RotationZ< T > ;
        using P  = qclab::qgates::Phase< T > ;
        using CX = qclab::qgates::CX< T > ;
        using CY = qclab::qgates::CY< T > ;
        using CZ = qclab::qgates::CZ< T > ;
        using SWAP = qclab::qgates::SWAP< T > ;

        // loop over commands
        std::unique_ptr< G > gate ;
        std::vector< std::unique_ptr< G > > gates ;
        try {
          for ( auto command : gateCommands_ ) {
            if ( parse1const< H >( "h" , command , gate ) ) {       // Hadamard
              gates.push_back( std::move( gate ) ) ;
            } else if ( parse1const< X >( "x" , command , gate ) ) {// PauliX
              gates.push_back( std::move( gate ) ) ;
            } else if ( parse1const< Y >( "y" , command , gate ) ) {// PauliY
              gates.push_back( std::move( gate ) ) ;
            } else if ( parse1const< Z >( "z" , command , gate ) ) {// PauliZ
              gates.push_back( std::move( gate ) ) ;
            } else if ( parse1a< RX >( "rx" , command , gate ) ) {  // RotationX
              gates.push_back( std::move( gate ) ) ;
            } else if ( parse1a< RY >( "ry" , command , gate ) ) {  // RotationY
              gates.push_back( std::move( gate ) ) ;
            } else if ( parse1a< RZ >( "rz" , command , gate ) ) {  // RotationZ
              gates.push_back( std::move( gate ) ) ;
            } else if ( parse1a< P >( "p" , command , gate ) ) {    // Phase
              gates.push_back( std::move( gate ) ) ;
            } else if ( parse2const< CX >( "cx" , command , gate ) ) { // CX
              gates.push_back( std::move( gate ) ) ;
            } else if ( parse2const< CY >( "cy" , command , gate ) ) { // CY
              gates.push_back( std::move( gate ) ) ;
            } else if ( parse2const< CZ >( "cz" , command , gate ) ) { // CZ
              gates.push_back( std::move( gate ) ) ;
            } else if ( parse2const< SWAP >( "swap" , command , gate ) ) {//SWAP
              gates.push_back( std::move( gate ) ) ;
            }
          }
        } catch ( const std::exception& e ) {
          std::cerr << "Failed parsing QASM file!" << std::endl ;
        }
        return gates ;
      }

      /// Parses a qubit.
      int parseQubit( std::string& str ) const {
        return read_value< int >( left_of( str , "]" ) ) ;
      }

      /// Parses an angle.
      template <typename T>
      T parseAngle( std::string& str ) const {
        return read_value< T >( left_of( str , ")" ) ) ;
      }

      /// Parses a constant 1-qubit gate.
      template <typename G, typename GATE>
      int parse1const( std::string name , std::string& command ,
                       std::unique_ptr< GATE >& gate ) const {
        const auto n1 = name.length() ;
        const auto n2 = qregName_.length() ;
        if ( command.substr( 0 , n1 + n2 ) == name + qregName_ ) {
          const auto qubit = parseQubit( command.erase( 0 , n1 + n2 + 1 ) ) ;
          gate = std::make_unique< G >( qubit ) ;
          return 1 ;
        }
        return 0 ;
      }

      /// Parses a 1-qubit gate with 1 parameter.
      template <typename G, typename GATE>
      int parse1a( std::string name , std::string& command ,
                   std::unique_ptr< GATE >& gate ) const {
        using T = typename G::value_type ;
        using R = qclab::real_t< T > ;
        const auto n1 = name.length() ;
        const auto n2 = qregName_.length() ;
        if ( command.substr( 0 , n1 + 1 ) == name + "(" ) {
          const auto angle = parseAngle< R >( command.erase( 0 , n1 + 1 ) ) ;
          const auto qubit = parseQubit(      command.erase( 0 , n2 + 1 ) ) ;
          gate = std::make_unique< G >( qubit , angle ) ;
          return 1 ;
        }
        return 0 ;
      }

      /// Parses a constant 2-qubit gate.
      template <typename G, typename GATE>
      int parse2const( std::string name , std::string& command ,
                       std::unique_ptr< GATE >& gate ) const {
        const auto n1 = name.length() ;
        const auto n2 = qregName_.length() ;
        if ( command.substr( 0 , n1 + n2 ) == name + qregName_ ) {
          const auto qubit0 = parseQubit( command.erase( 0 , n1 + n2 + 1 ) ) ;
          const auto qubit1 = parseQubit( command.erase( 0 , n2 + 2 ) ) ;
          gate = std::make_unique< G >( qubit0 , qubit1 ) ;
          return 1 ;
        }
        return 0 ;
      }

    private:
      /// Parses the QASM file.
      void parse() ;

      /// Input stream of this QASM file.
      std::ifstream               stream_ ;
      /// Number of qubits of this QASM file.
      int                         nbQubits_ ;
      /// Quantum register name of this QASM file.
      std::string                 qregName_ ;
      /// Gates of this QASM file.
      std::vector< std::string >  gateCommands_ ;

  } ; // class QASMFile

} // namespace qclab::io

