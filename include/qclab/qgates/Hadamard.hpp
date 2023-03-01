//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#pragma once

#include "qclab/qgates/QGate1.hpp"

namespace qclab {

  namespace qgates {

    /**
     * \class Hadamard
     * \brief 1-qubit Hadamard gate.
     */
    template <typename T>
    class Hadamard : public QGate1< T >
    {

      public:
        /// Default constructor. Constructs a Hadamard gate on qubit 0.
        Hadamard()
        : QGate1< T >()
        { } // Hadamard()

        /// Constructs a Hadamard gate on the given qubit `qubit`.
        Hadamard( const int qubit )
        : QGate1< T >( qubit )
        { } // Hadamard(qubit)

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
          using R = qclab::real_t< T > ;
          const R sqrt2 = R(1) / std::sqrt( R(2) ) ;
          return qclab::dense::SquareMatrix< T >( sqrt2 ,  sqrt2 ,
                                                  sqrt2 , -sqrt2 ) ;
        }

        // apply
        void apply( Op op , const int nbQubits , std::vector< T >& vector ,
                    const int offset = 0 ) const override ;

      #ifdef QCLAB_OMP_OFFLOADING
        // apply_device
        void apply_device( Op op , const int nbQubits , T* vector ,
                           const int offset = 0 ) const override ;
      #endif

        // apply
        void apply( Side side , Op op , const int nbQubits ,
                    qclab::dense::SquareMatrix< T >& matrix ,
                    const int offset = 0 ) const override ;

        // print

        // toQASM
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          stream << qasmH( this->qubit_ + offset ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using H = Hadamard< T > ;
          if ( const H* p = dynamic_cast< const H* >( &other ) ) return true ;
          return false ;
        }

    } ; // class Hadamard

  } // namespace qgates

} // namespace qclab

