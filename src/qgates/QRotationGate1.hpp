//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_QRotationGate1_hpp
#define qclab_qgates_QRotationGate1_hpp

#include "QGate1.hpp"
#include "QRotation.hpp"
#include <cassert>

namespace qclab {

  namespace qgates {

    /**
     * \class QRotationGate1
     * \brief Base class for 1-qubit rotation gates.
     */
    template <typename T>
    class QRotationGate1 : public QGate1< T > ,
                           public QRotation< qclab::real_t< T > >
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
        , QRotation< real_type >()
        { } // QRotationGate1()

        /**
         * \brief Constructs a 1-qubit rotation gate on the given qubit `qubit`
         *        with quantum angle `angle` = \f$\theta/2\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        QRotationGate1( const int qubit , const angle_type& angle ,
                        const bool fixed = false )
        : QGate1< T >( qubit )
        , QRotation< real_type >( angle , fixed )
        { } // QRotationGate1(qubit,angle,fixed)

        /**
         * \brief Constructs a 1-qubit rotation gate on the given qubit `qubit`
         *        with value `theta` = \f$\theta\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        QRotationGate1( const int qubit , const real_type theta ,
                        const bool fixed = false )
        : QGate1< T >( qubit )
        , QRotation< real_type >( theta , fixed )
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
        , QRotation< real_type >( cos , sin , fixed )
        { } // QRotationGate1(qubit,cos,sin,fixed)

        // nbQubits

        // fixed
        inline bool fixed() const override {
          return QRotation< real_type >::fixed() ;
        }

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

        // angle

        // theta

        // cos

        // sin

        // update(angle)

        // update(theta)

        // update(cos,sin)

    } ; // class QRotationGate1

  } // namespace qgates

} // namespace qclab

#endif

