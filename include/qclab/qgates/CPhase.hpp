//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_CPhase_hpp
#define qclab_qgates_CPhase_hpp

#include "qclab/qgates/QControlledGate2.hpp"
#include "qclab/qgates/Phase.hpp"

namespace qclab {

  namespace qgates {

    /**
     * \class CPhase
     * \brief Controlled phase gate.
     */
    template <typename T>
    class CPhase : public QControlledGate2< T >
    {

      public:
        /// Real value type of this controlled phase gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum angle type of this controlled phase gate.
        using angle_type = qclab::QAngle< real_type > ;
        /// Gate type of the controlled 1-qubit gate of this 2-qubit gate.
        using gate_type = Phase< T > ;

        /**
         * \brief Default constructor. Constructs a controlled phase gate with
         *        control qubit 0, target qubit 1, and parameter
         *        \f$\theta = 0\f$.
         */
        CPhase()
        : QControlledGate2< T >( 0 )
        , gate_( std::make_unique< Phase< T > >( 1 , 0. ) )
        { } // CPhase()

        /**
         * \brief Constructs a controlled phase gate with given control qubit
         *        `control`, target qubit `target`, quantum angle `angle` =
         *        \f$\theta\f$, and control state `controlState`.
         *        The default control state is 1.
         */
        CPhase( const int control , const int target , const angle_type& angle ,
                const int controlState = 1 )
        : QControlledGate2< T >( control , controlState )
        , gate_( std::make_unique< Phase< T > >( target , angle ) )
        {
          assert( control >= 0 ) ; assert( target >= 0 ) ;
          assert( control != target ) ;
        } // CPhase(control,target,angle,controlState)

        /**
         * \brief Constructs a controlled phase gate with given control qubit
         *        `control`, target qubit `target`, value `theta` =
         *        \f$\theta\f$, and control state `controlState`.
         *        The default control state is 1.
         */
        CPhase( const int control , const int target , const real_type theta ,
                const int controlState = 1 )
        : QControlledGate2< T >( control , controlState )
        , gate_( std::make_unique< Phase< T > >( target , theta ) )
        {
          assert( control >= 0 ) ; assert( target >= 0 ) ;
          assert( control != target ) ;
        } // CPhase(control,target,theta,controlState)

        /**
         * \brief Constructs a controlled phase gate with given control
         *        qubit `control`, target qubit `target`, values
         *        `cos` = \f$\cos(\theta)\f$ and `sin` = \f$\sin(\theta)\f$, and
         *        control state `controlState`. The default control state is 1.
         */
        CPhase( const int control , const int target ,
                const real_type cos , real_type sin ,
                const int controlState = 1 )
        : QControlledGate2< T >( control , controlState )
        , gate_( std::make_unique< Phase< T > >( target , cos , sin ) )
        {
          assert( control >= 0 ) ; assert( target >= 0 ) ;
          assert( control != target ) ;
        } // CPhase(control,target,cos,sin,controlState)

        // nbQubits

        // fixed
        inline bool fixed() const override { return gate_->fixed() ; }

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
          stream << qasmCPhase( this->control() + offset ,
                                this->target()  + offset , theta() ) ;
          if ( this->controlState() == 0 ) {
            stream << qasmX( this->control() + offset ) ;
          }
          return 0 ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using CP = CPhase< T > ;
          if ( const CP* p = dynamic_cast< const CP* >( &other ) ) {
            if ( p->angle() != this->angle() ) return false ;
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

        /// Makes this controlled phase gate fixed.
        inline void makeFixed() { gate_->makeFixed() ; }

        /// Makes this controlled phase gate variable.
        inline void makeVariable() { gate_->makeVariable() ; }

        /**
         * \brief Returns the quantum angle \f$\theta\f$ of this controlled
         *        phase gate.
         */
        inline const angle_type& angle() const {
          return gate_->angle() ;
        }

        /**
         * \brief Returns the numerical value \f$\theta\f$ of this controlled
         *        phase gate.
         */
        inline real_type theta() const {
          return gate_->theta() ;
        }

        /// Returns the cosine \f$\cos(\theta)\f$ of this controlled phase gate.
        inline real_type cos() const {
          return gate_->cos() ;
        }

        /// Returns the sine \f$\sin(\theta)\f$ of this controlled phase gate.
        inline real_type sin() const {
          return gate_->sin() ;
        }

        /**
         * \brief Updates this controlled phase gate with the given quantum
         *        angle `angle` = \f$\theta\f$.
         */
        void update( const angle_type& angle ) { gate_->update( angle ) ; }

        /**
         * \brief Updates this controlled phase gate with the given value
         *        `theta` = \f$\theta\f$.
         */
        void update( const real_type theta ) { gate_->update( theta ) ; }

        /**
         * \brief Updates this controlled phase gate with the given values
         *        `cos` = \f$\cos(\theta)\f$ and `sin` = \f$\sin(\theta)\f$.
         */
        void update( const real_type cos , const real_type sin ) {
          gate_->update( cos , sin ) ;
        }

      protected:
        /// Unique pointer to the phase gate of this controlled phase gate.
        std::unique_ptr< gate_type >  gate_ ;

    } ; // class CPhase

  } // namespace qgates

} // namespace qclab

#endif

