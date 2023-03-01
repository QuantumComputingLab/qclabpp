//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#pragma once

#include "qclab/qgates/QGate1.hpp"

namespace qclab {

  namespace qgates {

    /**
     * \class Phase90
     * \brief 1-qubit S gate (pi/2 rotation about the Z-axis).
     */
    template <typename T>
    class Phase90 : public QGate1< T >
    {

      public:
        /// Default constructor. Constructs a S gate on qubit 0.
        Phase90()
        : QGate1< T >()
        { } // Phase90()

        /// Constructs a S gate on the given qubit `qubit`.
        Phase90( const int qubit )
        : QGate1< T >( qubit )
        { } // Phase90(qubit)

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
          return qclab::dense::SquareMatrix< T >( 1 ,   0    ,
                                                  0 , T(0,1) ) ;
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
          stream << qasmS( this->qubit_ + offset ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using S = Phase90< T > ;
          if ( const S* p = dynamic_cast< const S* >( &other ) ) return true ;
          return false ;
        }

    } ; // class Phase90

  } // namespace qgates

} // namespace qclab

