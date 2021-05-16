//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_PauliX_hpp
#define qclab_qgates_PauliX_hpp

#include "QGate1.hpp"
#include <cassert>

namespace qclab {

  namespace qgates {

    /**
     * \class PauliX
     * \brief 1-qubit Pauli-X gate.
     */
    template <typename T>
    class PauliX : public QGate1< T >
    {

      public:
        /// Default constructor. Constructs a Pauli-X gate on qubit 0.
        PauliX()
        : QGate1< T >()
        { } // PauliX()

        /// Constructs a Pauli-X gate on the given qubit `qubit`.
        PauliX( const int qubit )
        : QGate1< T >( qubit )
        { } // PauliX(qubit)

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
          return qclab::dense::SquareMatrix< T >( 0 , 1 ,
                                                  1 , 0 ) ;
        }

        // apply

        // print

        // toQASM
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          stream << qasmX( this->qubit_ + offset ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using X = PauliX< T > ;
          if ( const X* p = dynamic_cast< const X* >( &other ) ) return true ;
          return false ;
        }

    } ; // class PauliX

  } // namespace qgates

} // namespace qclab

#endif

