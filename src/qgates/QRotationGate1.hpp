//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_QRotationGate1_hpp
#define qclab_qgates_QRotationGate1_hpp

#include "QGate1.hpp"
#include "QAdjustable.hpp"
#include "QRotation.hpp"
#include <cassert>

namespace qclab {

  namespace qgates {

    /**
     * \class QRotationGate1
     * \brief Base class for 1-qubit rotation gates.
     */
    template <typename T>
    class QRotationGate1 : public QGate1< T > , public QAdjustable
    {

      public:
        /// Real value type of this 1-qubit rotation gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum angle type of this 1-qubit rotation gate.
        using angle_type = qclab::QAngle< real_type > ;

        /**
         * \brief Default constructor. Constructs a 1-qubit rotation gate on
         *        qubit 0 with parameter \f$\theta = 0\f$.
         */
        QRotationGate1()
        : QGate1< T >( 0 )
        , rotation_()
        , QAdjustable()
        { } // QRotationGate1()

        /**
         * \brief Constructs a 1-qubit rotation gate on the given qubit `qubit`
         *        with quantum angle `angle` = \f$\theta/2\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        QRotationGate1( const int qubit , const angle_type& angle ,
                        const bool fixed = false )
        : QGate1< T >( qubit )
        , rotation_( angle )
        , QAdjustable( fixed )
        { } // QRotationGate1(qubit,angle,fixed)

        /**
         * \brief Constructs a 1-qubit rotation gate on the given qubit `qubit`
         *        with value `theta` = \f$\theta\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        QRotationGate1( const int qubit , const real_type theta ,
                        const bool fixed = false )
        : QGate1< T >( qubit )
        , rotation_( theta )
        , QAdjustable( fixed )
        { } // QRotationGate1(qubit,theta,fixed)

        /**
         * \brief Constructs a 1-qubit rotation gate on the given qubit `qubit`
         *        with values `cos` = \f$\cos(\theta/2)\f$ and
         *        `sin` = \f$\sin(\theta/2)\f$, and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        QRotationGate1( const int qubit , const real_type cos ,
                        const real_type sin , const bool fixed = false )
        : QGate1< T >( qubit )
        , rotation_( cos , sin )
        , QAdjustable( fixed )
        { } // QRotationGate1(qubit,cos,sin,fixed)

        // nbQubits

        // fixed
        inline bool fixed() const override { return QAdjustable::fixed() ; }

        // controlled

        // qubit

        // setQubit

        // qubits

        // setQubits

        // matrix

        // apply

        // print

        // toQASM

        /// Checks if `other` is equal to this 1-qubit rotation gate.
        inline bool operator==( const QRotationGate1< T >& other ) const {
          return other.angle() == this->angle() ;
        }

        /// Checks if `other` is different from this 1-qubit rotation gate.
        inline bool operator!=( const QRotationGate1< T >& other ) const {
          return !( *this == other ) ;
        }

        // equals

        /// Returns the quantum angle \f$\theta/2\f$ of this rotation gate.
        inline const angle_type& angle() const { return rotation_.angle() ; }

        /// Returns the numerical value \f$\theta\f$ of this rotation gate.
        inline real_type theta() const { return rotation_.theta() ; }

        /// Returns the cosine \f$\cos(\theta/2)\f$ of this rotation gate.
        inline real_type cos() const { return rotation_.cos() ; }

        /// Returns the sine \f$\sin(\theta/2)\f$ of this rotation gate.
        inline real_type sin() const { return rotation_.sin() ; }

        /**
         * \brief Updates this rotation gate with the given quantum angle
         *        `angle` = \f$\theta/2\f$.
         */
        void update( const angle_type& angle ) {
          assert( !fixed() ) ;
          rotation_.update( angle ) ;
        }

        /**
         * \brief Updates this rotation gate with the given value
         *        `theta` = \f$\theta\f$.
         */
        void update( const real_type theta ) {
          assert( !fixed() ) ;
          rotation_.update( theta ) ;
        }

        /**
         * \brief Updates this rotation gate with the given values
         *        `cos` = \f$\cos(\theta/2)\f$ and `sin` = \f$\sin(\theta/2)\f$.
         */
        void update( const real_type cos , const real_type sin ) {
          assert( !fixed() ) ;
          rotation_.update( cos , sin ) ;
        }

      protected:
        /// Quantum rotation of this rotation gate.
        QRotation< real_type >  rotation_ ;

    } ; // class QRotationGate1

  } // namespace qgates

} // namespace qclab

#endif

