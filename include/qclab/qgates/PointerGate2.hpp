//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#pragma once

#include "qclab/qgates/QGate2.hpp"

namespace qclab {

  namespace qgates {

    /**
     * \class PointerGate2
     * \brief 2-qubit gate pointing to another 2-qubit gate.
     */
    template <typename T>
    class PointerGate2 : public QGate2< T >
    {

      public:
        /**
         * \brief Constructs a 2-qubit pointer gate from the given gate pointer
         *        `gate` and qubit offset `offset`.
         *        The default value of `offset` is 0.
         */
        PointerGate2( QGate2< T >* gate , const int offset = 0 )
        : offset_( offset )
        , gate_( gate )
        { } // PointerGate2(gate,offset)

        // nbQubits
        inline int nbQubits() const override { return 2 ; }

        // fixed
        inline bool fixed() const override { return gate_->fixed() ; }

        // controlled
        inline bool controlled() const override {
          return gate_->controlled() ;
        }

        // qubit
        inline int qubit() const override { return gate_->qubit() + offset_ ; }

        // setQubit

        // qubits
        std::vector< int > qubits() const override {
          auto v = gate_->qubits() ;
          for ( size_t i = 0; i < v.size(); ++i ) { v[i] += offset_ ; }
          return v ;
        }

        // setQubits
        inline void setQubits( const int* qubits ) override { assert( false ) ;}

        // matrix
        qclab::dense::SquareMatrix< T > matrix() const override {
          return gate_->matrix() ;
        }

        // apply
        void apply( Op op , const int nbQubits , std::vector< T >& vector ,
                    const int offset = 0 ) const override {
          gate_->apply( op , nbQubits , vector , offset_ + offset ) ;
        }

      #ifdef QCLAB_OMP_OFFLOADING
        // apply_device
        void apply_device( Op op , const int nbQubits , T* vector ,
                           const int offset = 0 ) const override {
          gate_->apply_device( op , nbQubits , vector , offset_ + offset ) ;
        }
      #endif

        // apply
        void apply( Side side , Op op , const int nbQubits ,
                    qclab::dense::SquareMatrix< T >& matrix ,
                    const int offset = 0 ) const override {
          gate_->apply( side , op , nbQubits , matrix , offset_ + offset ) ;
        }

        // print

        // toQASM
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          return gate_->toQASM( stream , offset_ + offset ) ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using G = PointerGate2< T > ;
          if ( const G* p = dynamic_cast< const G* >( &other ) ) {
            return *p->ptr() == *this->ptr() ;
          }
          return other == *this->ptr() ;
        }

        /// Returns the qubit offset of this 2-qubit pointer gate.
        int offset() const { return offset_ ; }

        /// Sets the qubit offset of this 2-qubit pointer gate.
        void setOffset( const int offset ) { offset_ = offset ; }

        /// Returns a pointer to the gate of this 2-qubit pointer gate.
        QGate2< T >* ptr() { return gate_ ; }

        /// Returns a const pointer to the gate of this 2-qubit pointer gate.
        const QGate2< T >* ptr() const { return gate_ ; }

      protected:
        int           offset_ ;  ///< Qubit offset of this 2-qubit pointer gate.
        QGate2< T >*  gate_ ;    ///< Gate pointer of this 2-qubit pointer gate.

    } ; // class PointerGate2

  } // namespace qgates

} // namespace qclab

