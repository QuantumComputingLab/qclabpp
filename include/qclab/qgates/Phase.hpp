//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_Phase_hpp
#define qclab_qgates_Phase_hpp

#include "qclab/qgates/QGate1.hpp"
#include "qclab/QAdjustable.hpp"
#include "qclab/QAngle.hpp"

namespace qclab {

  namespace qgates {

    /**
     * \class Phase
     * \brief 1-qubit phase gate.
     */
    template <typename T>
    class Phase : public QGate1< T > , public QAdjustable
    {

      public:
        /// Real value type of this phase gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum angle type of this phase gate.
        using angle_type = qclab::QAngle< real_type > ;

        /**
         * \brief Default constructor. Constructs a phase gate on qubit 0 with
         *        parameter \f$\theta = 0\f$.
         */
        Phase()
        : QGate1< T >( 0 )
        , angle_( 0 )
        , QAdjustable( false )
        { } // Phase()

        /**
         * \brief Constructs a phase gate on qubit 0 with the given quantum
         *        angle `angle` = \f$\theta\f$.
         */
        Phase( const angle_type& angle )
        : QGate1< T >( 0 )
        , angle_( angle )
        , QAdjustable( false )
        { } // Phase(angle)

        /**
         * \brief Constructs a phase gate on qubit 0 with the given value
         *        `theta` = \f$\theta\f$.
         */
        Phase( const real_type theta )
        : QGate1< T >( 0 )
        , angle_( theta )
        , QAdjustable( false )
        { } // Phase(theta)

        /**
         * \brief Constructs a phase gate on qubit 0 with the given values
         *        `cos` = \f$\cos(\theta)\f$ and `sin` = \f$\sin(\theta)\f$.
         */
        Phase( const real_type cos , const real_type sin )
        : QGate1< T >( 0 )
        , angle_( cos , sin )
        , QAdjustable( false )
        { } // Phase(cos,sin)

        /**
         * \brief Constructs a phase gate on the given qubit `qubit` with
         *        quantum angle `angle` = \f$\theta\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        Phase( const int qubit , const angle_type& angle ,
               const bool fixed = false )
        : QGate1< T >( qubit )
        , angle_( angle )
        , QAdjustable( fixed )
        { } // Phase(qubit,angle,fixed)

        /**
         * \brief Constructs a phase gate on the given qubit `qubit` with
         *        value `theta` = \f$\theta\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        Phase( const int qubit , const real_type theta ,
               const bool fixed = false )
        : QGate1< T >( qubit )
        , angle_( theta )
        , QAdjustable( fixed )
        { } // Phase(qubit,theta,fixed)

        /**
         * \brief Constructs a phase gate on the given qubit `qubit` with values
         *        `cos` = \f$\cos(\theta)\f$ and `sin` = \f$\sin(\theta)\f$, and
         *        flag `fixed`. The default value of `fixed` is false.
         */
        Phase( const int qubit , const real_type cos , const real_type sin ,
               const bool fixed = false )
        : QGate1< T >( qubit )
        , angle_( cos , sin )
        , QAdjustable( fixed )
        { } // Phase(qubit,cos,sin,fixed)

        // nbQubits

        // fixed
        inline bool fixed() const override { return QAdjustable::fixed() ; }

        // controlled

        // qubit

        // setQubit

        // qubits

        // setQubits

        // matrix
        qclab::dense::SquareMatrix< T > matrix() const override {
          return qclab::dense::SquareMatrix< T >( T(1) , T(0) ,
                                                  T(0) , T( cos() , sin() ) ) ;
        }

        // apply

        // print

        // toQASM
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          stream << qasmRz( this->qubit_ + offset , theta() ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using Z = Phase< T > ;
          if ( const Z* p = dynamic_cast< const Z* >( &other ) ) {
            return ( p->angle() == this->angle() ) ;
          }
          return false ;
        }

        /// Returns the quantum angle \f$\theta\f$ of this phase gate.
        inline const angle_type& angle() const { return angle_ ; }

        /// Returns the numerical value \f$\theta\f$ of this phase gate.
        inline real_type theta() const { return angle_.theta() ; }

        /// Returns the cosine \f$\cos(\theta)\f$ of this phase gate.
        inline real_type cos() const { return angle_.cos() ; }

        /// Returns the sine \f$\sin(\theta)\f$ of this phase gate.
        inline real_type sin() const { return angle_.sin() ; }

        /**
         * \brief Updates this phase gate with the given quantum angle
         *        `angle` = \f$\theta\f$.
         */
        void update( const angle_type& angle ) {
          assert( !fixed() ) ;
          angle_.update( angle ) ;
        }

        /// Updates this phase gate with the given value `theta` = \f$\theta\f$.
        void update( const real_type theta ) {
          assert( !fixed() ) ;
          angle_.update( theta ) ;
        }

        /**
         * \brief Updates this phase gate with the given values
         *        `cos` = \f$\cos(\theta)\f$ and `sin` = \f$\sin(\theta)\f$.
         */
        void update( const real_type cos , const real_type sin ) {
          assert( !fixed() ) ;
          angle_.update( cos , sin ) ;
        }

      protected:
        angle_type  angle_ ;  ///< Quantum angle of this phase gate.

    } ; // class Phase

  } // namespace qgates

} // namespace qclab

#endif

