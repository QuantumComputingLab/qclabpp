//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_PauliY_hpp
#define qclab_qgates_PauliY_hpp

#include "QGate1.hpp"
#include <cassert>

namespace qclab {

  namespace qgates {

    /**
     * \class PauliY
     * \brief 1-qubit Pauli-Y gate.
     */
    template <typename T>
    class PauliY : public QGate1< T >
    {

      public:
        /// Default constructor. Constructs a Pauli-Y gate on qubit 0.
        PauliY()
        : QGate1< T >()
        { } // PauliY()

        /// Constructs a Pauli-Y gate on the given qubit `qubit`.
        PauliY( const int qubit )
        : QGate1< T >( qubit )
        { } // PauliY(qubit)

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
          return qclab::dense::SquareMatrix< T >(   0    , T(0,-1) ,
                                                  T(0,1) ,   0     ) ;
        }

        // apply

        // print

        // toQASM
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          stream << qasmY( this->qubit_ + offset ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using Y = PauliY< T > ;
          if ( const Y* p = dynamic_cast< const Y* >( &other ) ) return true ;
          return false ;
        }

    } ; // class PauliY

  } // namespace qgates

} // namespace qclab

#endif

