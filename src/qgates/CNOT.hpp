//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_CNOT_hpp
#define qclab_qgates_CNOT_hpp

#include "QControlledGate2.hpp"
#include "PauliX.hpp"
#include <cassert>
#include <vector>

namespace qclab {

  namespace qgates {

    /**
     * \class CNOT
     * \brief Controlled NOT gate (CNOT).
     */
    template <typename T>
    class CNOT : public QControlledGate2< T >
    {

      public:
        /// Gate type of the controlled 1-qubit gate of this 2-qubit gate.
        using gate_type = PauliX< T > ;

        /**
         * \brief Default constructor. Constructs a CNOT gate with
         *        control qubit 0 and target qubit 1.
         */
        CNOT()
        : QControlledGate2< T >( 0 )
        , gate_( std::make_unique< PauliX< T > >( 1 ) )
        { } // CNOT()

        /**
         * \brief Constructs a CNOT gate with given control qubit `control`,
         *        target qubit `target`, and control state `controlState`.
         *        The default control state is 1.
         */
        CNOT( const int control , const int target ,
              const int controlState = 1 )
        : QControlledGate2< T >( control , controlState )
        , gate_( std::make_unique< PauliX< T > >( target ) )
        {
          assert( control >= 0 ) ; assert( target >= 0 ) ;
          assert( control != target ) ;
        } // CNOT(control,target,controlState)

        // nbQubits

        // fixed
        inline bool fixed() const override { return true ; }

        // controlled

        // qubit
        inline int qubit() const override {
          return std::min( this->control_ , gate_->qubit() ) ;
        }

        // setQubit

        // qubits

        // setQubits

        // matrix

        // apply

        // print

        // toQASM
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          if ( this->controlState() == 0 ) {
            stream << qasmX( this->control() + offset ) ;
          }
          stream << qasmCNOT( this->control() + offset ,
                              this->target()  + offset ) ;
          if ( this->controlState() == 0 ) {
            stream << qasmX( this->control() + offset ) ;
          }
          return 0 ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using CX = CNOT< T > ;
          if ( const CX* p = dynamic_cast< const CX* >( &other ) ) {
            if ( p->controlState() != this->controlState() ) return false ;
            return ( ( p->control() < p->target() ) &&
                     ( this->control() < this->target() ) ) ||
                   ( ( p->control() > p->target() ) &&
                     ( this->control() > this->target() ) ) ;
          }
          return false ;
        }

        // control

        // target
        inline int target() const override { return gate_->qubit() ; }

        // controlState

        // gate
        inline const gate_type* gate() const override { return gate_.get() ; }

        // setControl

        // setTarget
        void setTarget( const int target ) override {
          assert( target >= 0 ) ; assert( this->control() != target ) ;
          gate_->setQubit( target ) ;
        }

        // setControlState

      protected:
        /// Unique pointer to the Pauli-X gate of this CNOT gate.
        std::unique_ptr< gate_type >  gate_ ;

    } ; // class CNOT

  } // namespace qgates

} // namespace qclab

#endif

