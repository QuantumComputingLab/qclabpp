//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_CRotationZ_hpp
#define qclab_qgates_CRotationZ_hpp

#include "qclab/qgates/QControlledGate2.hpp"
#include "qclab/qgates/RotationZ.hpp"

namespace qclab {

  namespace qgates {

    /**
     * \class CRotationZ
     * \brief Controlled Z-rotation gate.
     */
    template <typename T>
    class CRotationZ : public QControlledGate2< T >
    {

      public:
        /// Real value type of this controlled Z-rotation gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum rotation type of this controlled Z-rotation gate.
        using rotation_type = qclab::QRotation< real_type > ;
        /// Gate type of the controlled 1-qubit gate of this 2-qubit gate.
        using gate_type = RotationZ< T > ;

        /**
         * \brief Default constructor. Constructs a controlled Z-rotation gate
         *        with control qubit 0, target qubit 1, and parameter
         *        \f$\theta = 0\f$.
         */
        CRotationZ()
        : QControlledGate2< T >( 0 )
        , gate_( std::make_unique< RotationZ< T > >( 1 , 0. ) )
        { } // CRotationZ()

        /**
         * \brief Constructs a controlled Z-rotation gate with given control
         *        qubit `control`, target qubit `target`, quantum rotation
         *        `rot` = \f$\theta\f$, and control state `controlState`.
         *        The default control state is 1.
         */
        CRotationZ( const int control , const int target ,
                    const rotation_type& rot , const int controlState = 1 )
        : QControlledGate2< T >( control , controlState )
        , gate_( std::make_unique< RotationZ< T > >( target , rot ) )
        {
          assert( control >= 0 ) ; assert( target >= 0 ) ;
          assert( control != target ) ;
        } // CRotationZ(control,target,rot,controlState)

        /**
         * \brief Constructs a controlled Z-rotation gate with given control
         *        qubit `control`, target qubit `target`, value `theta` =
         *        \f$\theta\f$, and control state `controlState`.
         *        The default control state is 1.
         */
        CRotationZ( const int control , const int target ,
                    const real_type theta , const int controlState = 1 )
        : QControlledGate2< T >( control , controlState )
        , gate_( std::make_unique< RotationZ< T > >( target , theta ) )
        {
          assert( control >= 0 ) ; assert( target >= 0 ) ;
          assert( control != target ) ;
        } // CRotationZ(control,target,theta,controlState)

        /**
         * \brief Constructs a controlled Z-rotation gate with given control
         *        qubit `control`, target qubit `target`, values
         *        `cos` = \f$\cos(\theta/2)\f$ and `sin` = \f$\sin(\theta/2)\f$,
         *        and control state `controlState`.
         *        The default control state is 1.
         */
        CRotationZ( const int control , const int target ,
                    const real_type cos , real_type sin ,
                    const int controlState = 1 )
        : QControlledGate2< T >( control , controlState )
        , gate_( std::make_unique< RotationZ< T > >( target , cos , sin ) )
        {
          assert( control >= 0 ) ; assert( target >= 0 ) ;
          assert( control != target ) ;
        } // CRotationZ(control,target,cos,sin,controlState)

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
          stream << qasmCRz( this->control() + offset ,
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
          using CRZ = CRotationZ< T > ;
          if ( const CRZ* p = dynamic_cast< const CRZ* >( &other ) ) {
            if ( p->rotation() != this->rotation() ) return false ;
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

        /// Makes this controlled Z-rotation gate fixed.
        inline void makeFixed() { gate_->makeFixed() ; }

        /// Makes this controlled Z-rotation gate variable.
        inline void makeVariable() { gate_->makeVariable() ; }

        /**
         * \brief Returns the quantum rotation \f$\theta\f$ of this controlled
         *        Z-rotation gate.
         */
        inline const rotation_type& rotation() const {
          return gate_->rotation() ;
        }

        /**
         * \brief Returns the numerical value \f$\theta\f$ of this controlled
         *        Z-rotation gate.
         */
        inline real_type theta() const {
          return gate_->theta() ;
        }

        /**
         * \brief Returns the cosine \f$\cos(\theta/2)\f$ of this controlled
         *        Z-rotation gate.
         */
        inline real_type cos() const {
          return gate_->cos() ;
        }

        /**
         * \brief Returns the sine \f$\sin(\theta/2)\f$ of this controlled
         *        Z-rotation gate.
         */
        inline real_type sin() const {
          return gate_->sin() ;
        }

        /**
         * \brief Updates this controlled Z-rotation gate with the given quantum
         *        rotation `rot` = \f$\theta\f$.
         */
        void update( const rotation_type& rot ) { gate_->update( rot ) ; }

        /**
         * \brief Updates this controlled Z-rotation gate with the given value
         *        `theta` = \f$\theta\f$.
         */
        void update( const real_type theta ) { gate_->update( theta ) ; }

        /**
         * \brief Updates this controlled Z-rotation gate with the given values
         *        `cos` = \f$\cos(\theta/2)\f$ and `sin` = \f$\sin(\theta/2)\f$.
         */
        void update( const real_type cos , const real_type sin ) {
          gate_->update( cos , sin ) ;
        }

      protected:
        /**
         * \brief Unique pointer to the Z-rotation gate of this controlled
         *        Z-rotation gate.
         */
        std::unique_ptr< gate_type >  gate_ ;

    } ; // class CRotationZ

  } // namespace qgates

} // namespace qclab

#endif

