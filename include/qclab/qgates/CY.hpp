//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_CY_hpp
#define qclab_qgates_CY_hpp

#include "qclab/qgates/QControlledGate2.hpp"
#include "qclab/qgates/PauliY.hpp"

namespace qclab {

  namespace qgates {

    /**
     * \class CY
     * \brief Controlled Pauli-Y gate.
     */
    template <typename T>
    class CY : public QControlledGate2< T >
    {

      public:
        /// Gate type of the controlled 1-qubit gate of this 2-qubit gate.
        using gate_type = PauliY< T > ;

        /**
         * \brief Default constructor. Constructs a controlled Pauli-Y gate
         *        with control qubit 0 and target qubit 1.
         */
        CY()
        : QControlledGate2< T >( 0 )
        , gate_( std::make_unique< PauliY< T > >( 1 ) )
        { } // CY()

        /**
         * \brief Constructs a controlled Pauli-Y gate with given control qubit
         *        `control`, target qubit `target`, and control state
         *        `controlState`. The default control state is 1.
         */
        CY( const int control , const int target ,
                 const int controlState = 1 )
        : QControlledGate2< T >( control , controlState )
        , gate_( std::make_unique< PauliY< T > >( target ) )
        {
          assert( control >= 0 ) ; assert( target >= 0 ) ;
          assert( control != target ) ;
        } // CY(control,target,controlState)

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
          stream << qasmCY( this->control() + offset ,
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
          using CG = CY< T > ;
          if ( const CG* p = dynamic_cast< const CG* >( &other ) ) {
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
        /// Unique pointer to the Pauli-Y gate of this CY gate.
        std::unique_ptr< gate_type >  gate_ ;

    } ; // class CY

  } // namespace qgates

} // namespace qclab

#endif

