//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#pragma once

#include "qclab/qgates/QGate2.hpp"
#include <array>

namespace qclab {

  namespace qgates {

    /**
     * \class SWAP
     * \brief SWAP gate.
     */
    template <typename T>
    class SWAP : public QGate2< T >
    {

      public:
        /// Default constructor. Constructs a SWAP gate between qubits 0 and 1.
        SWAP()
        : qubits_( { 0 , 1 } )
        {} // SWAP()

        /// Constructs a SWAP gate that swaps the qubits `qubit0` and `qubit1`.
        SWAP( const int qubit0 , const int qubit1 )
        : qubits_( { std::min( qubit0 , qubit1 ) ,
                     std::max( qubit0 , qubit1 ) } )
        {
          const int qubits[] = { qubit0 , qubit1 } ;
          setQubits( &qubits[0] ) ;
        } // SWAP(qubit0,qubit1)

        /// Constructs a SWAP gate that swaps the given qubits `qubits`.
        SWAP( const int* qubits )
        {
          setQubits( qubits ) ;
        } // SWAP(qubits)

        // nbQubits

        // fixed
        inline bool fixed() const override { return true ; }

        // controlled
        inline bool controlled() const override { return false ; }

        // qubit
        inline int qubit() const override { return qubits_[0] ; }

        // setQubit

        // qubits
        std::vector< int > qubits() const override {
          return std::vector< int >( { qubits_[0] , qubits_[1] } ) ;
        }

        // setQubits
        inline void setQubits( const int* qubits ) override {
          assert( qubits[0] >= 0 ) ; assert( qubits[1] >= 0 ) ;
          assert( qubits[0] != qubits[1] ) ;
          qubits_[0] = std::min( qubits[0] , qubits[1] ) ;
          qubits_[1] = std::max( qubits[0] , qubits[1] ) ;
        }

        // matrix
        qclab::dense::SquareMatrix< T > matrix() const override {
          return qclab::dense::SquareMatrix< T >( 1 , 0 , 0 , 0 ,
                                                  0 , 0 , 1 , 0 ,
                                                  0 , 1 , 0 , 0 ,
                                                  0 , 0 , 0 , 1 ) ;
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
          stream << qasmSWAP( qubits_[0] + offset , qubits_[1] + offset ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using G = SWAP< T > ;
          if ( const G* p = dynamic_cast< const G* >( &other ) ) return true ;
          return false ;
        }

        // control

        // target

        // controlState

        // gate

        // setControl

        // setTarget

        // setControlState

      protected:
        std::array< int , 2 >  qubits_ ;  ///< Qubits of this SWAP gate.

    } ; // class SWAP

  } // namespace qgates

} // namespace qclab

