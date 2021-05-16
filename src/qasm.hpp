//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qasm_hpp
#define qclab_qasm_hpp

#include <cassert>
#include <iostream>

namespace qclab {

  /// Returns a qasm sting for the given angle `angle`.
  inline auto qasm( const float angle ) {
    char buffer[16] ;
    int size = std::sprintf( buffer , "%.7f" , angle ) ;
    return std::string( buffer , size ) ;
  }

  /// Returns a qasm sting for the given angle `angle`.
  inline auto qasm( const double angle ) {
    char buffer[32] ;
    int size = std::sprintf( buffer , "%.15f" , angle ) ;
    return std::string( buffer , size ) ;
  }

  /// Returns a qasm string for a 1 qubit gate.
  inline auto qasm1( const char type[] , const int qubit ) {
    std::stringstream stream ;
    stream << type << " q[" << qubit << "];\n" ;
    return stream.str() ;
  }

  /// Returns a qasm string for a 2 qubit gate.
  inline auto qasm2( const char type[] , const int qubit0 , const int qubit1 ) {
    std::stringstream stream ;
    stream << type << " q[" << qubit0 << "], q[" << qubit1 << "];\n" ;
    return stream.str() ;
  }

  /// Returns a qasm string for a 1 qubit gate with angle `angle`.
  template <typename T>
  inline auto qasm1( const char type[] , const int qubit , const T angle ) {
    std::stringstream stream ;
    stream << type << "(" << qasm( angle ) << ") q[" << qubit << "];\n" ;
    return stream.str() ;
  }

  /// Returns a qasm string for a 2 qubit gate with angle `angle`.
  template <typename T>
  inline auto qasm2( const char type[] , const int qubit0 , const int qubit1 ,
                     const T angle ) {
    std::stringstream stream ;
    stream << type << "(" << qasm( angle ) << ") q[" << qubit0 << "], q["
                                                     << qubit1 << "];\n" ;
    return stream.str() ;
  }


  /// Returns a qasm string for an identity gate on the given qubit `qubit`.
  inline auto qasmI( const int qubit ) { return qasm1( "i" , qubit ) ; }

  /// Returns a qasm string for a Hadamard gate on the given qubit `qubit`.
  inline auto qasmH( const int qubit ) { return qasm1( "h" , qubit ) ; }

  /// Returns a qasm string for a Pauli-X gate on the given qubit `qubit`.
  inline auto qasmX( const int qubit ) { return qasm1( "x" , qubit ) ; }

  /// Returns a qasm string for a Pauli-Y gate on the given qubit `qubit`.
  inline auto qasmY( const int qubit ) { return qasm1( "y" , qubit ) ; }

  /// Returns a qasm string for a Pauli-Z gate on the given qubit `qubit`.
  inline auto qasmZ( const int qubit ) { return qasm1( "z" , qubit ) ; }


  /**
   * \brief Returns a qasm string for a Rotation X gate on the given qubit
   *        `qubit` and with angle `angle`.
   */
  template <typename T>
  inline auto qasmRx( const int qubit , const T angle ) {
    return qasm1( "rx" , qubit , angle ) ;
  }

  /**
   * \brief Returns a qasm string for a Rotation Y gate on the given qubit
   *        `qubit` and with angle `angle`.
   */
  template <typename T>
  inline auto qasmRy( const int qubit , const T angle ) {
    return qasm1( "ry" , qubit , angle ) ;
  }

  /**
   * \brief Returns a qasm string for a Rotation Z gate on the given qubit
   *        `qubit` and with angle `angle`.
   */
  template <typename T>
  inline auto qasmRz( const int qubit , const T angle ) {
    return qasm1( "rz" , qubit , angle ) ;
  }


  /**
   * \brief Returns a qasm string for a Rotation XX gate on the given qubits
   *        `qubit0` and `qubit1`, and with angle `angle`.
   */
  template <typename T>
  inline auto qasmRxx( const int qubit0 , const int qubit1 , const T angle ) {
    return qasm2( "rxx" , qubit0 , qubit1 , angle ) ;
  }

  /**
   * \brief Returns a qasm string for a Rotation YY gate on the given qubits
   *        `qubit0` and `qubit1`, and with angle `angle`.
   */
  template <typename T>
  inline auto qasmRyy( const int qubit0 , const int qubit1 , const T angle ) {
    return qasm2( "ryy" , qubit0 , qubit1 , angle ) ;
  }

  /**
   * \brief Returns a qasm string for a Rotation ZZ gate on the given qubits
   *        `qubit0` and `qubit1`, and with angle `angle`.
   */
  template <typename T>
  inline auto qasmRzz( const int qubit0 , const int qubit1 , const T angle ) {
    return qasm2( "rzz" , qubit0 , qubit1 , angle ) ;
  }


  /**
   * \brief Returns a qasm string for a CNOT gate with given control qubit
   *        `control` and target qubit `target`.
   */
  inline auto qasmCNOT( const int control , const int target ) {
    return qasm2( "cx" , control , target ) ;
  }

  /**
   * \brief Returns a qasm string for a controlled phase gate with given control
   *        qubit `control`, target qubit `target`, and angle `angle`.
   */
  template <typename T>
  inline auto qasmCPhase( const int control , const int target ,
                          const T angle ) {
    return qasm2( "cp" , control , target , angle ) ;
  }

  /**
   * \brief Returns a qasm string for a SWAP gate on the given qubits `qubit0`
   *        and `qubit1`.
   */
  inline auto qasmSWAP( const int qubit0 , const int qubit1 ) {
    return qasm2( "swap" , qubit0 , qubit1 ) ;
  }

} // namespace qclab

#endif

