//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_PauliZ_hpp
#define qclab_qgates_PauliZ_hpp

#include "QGate1.hpp"
#include <cassert>

namespace qclab {

  namespace qgates {

    /**
     * \class PauliZ
     * \brief 1-qubit Pauli-Z gate.
     */
    template <typename T>
    class PauliZ : public QGate1< T >
    {

      public:
        /// Default constructor. Constructs a Pauli-Z gate on qubit 0.
        PauliZ()
        : QGate1< T >()
        { } // PauliZ()

        /// Constructs a Pauli-Z gate on the given qubit `qubit`.
        PauliZ( const int qubit )
        : QGate1< T >( qubit )
        { } // PauliZ(qubit)

        // nbQubits

        // fixed
        inline bool fixed() const override { return true ; } ;

        // controlled

        // qubit

        // setQubit

        // qubits

        // setQubits

        // matrix
        qclab::dense::SquareMatrix< T > matrix() const override {
          return qclab::dense::SquareMatrix< T >( 1 ,  0 ,
                                                  0 , -1 ) ;
        }

        // apply

        // print

        // toQASM
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          stream << qasmZ( this->qubit_ + offset ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using Z = PauliZ< T > ;
          if ( const Z* p = dynamic_cast< const Z* >( &other ) ) return true ;
          return false ;
        }

    } ; // class PauliZ

  } // namespace qgates

} // namespace qclab

#endif

