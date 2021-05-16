//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_QRotationGate2_hpp
#define qclab_qgates_QRotationGate2_hpp

#include "QGate2.hpp"
#include "QRotation.hpp"
#include <cassert>

namespace qclab {

  namespace qgates {

    /**
     * \class QRotationGate2
     * \brief Base class for 2-qubit rotation gates.
     */
    template <typename T>
    class QRotationGate2 : public QGate2< T > ,
                           public QRotation< qclab::real_t< T > >
    {

      public:
        /// Real value type of this 2-qubit rotation gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum angle type of this 2-qubit rotation gate.
        using angle_type = qclab::QAngle< real_type > ;

        /**
         * \brief Default constructor. Constructs a 2-qubit rotation gate on
         *        qubits 0 and 1 with parameter \f$\theta = 0\f$.
         */
        QRotationGate2()
        : qubits_( { 0 , 1 } )
        , QRotation< real_type >()
        { } // QRotationGate2()

        /**
         * \brief Constructs a 2-qubit rotation gate on the given qubits
         *        `qubits` with quantum angle `angle` = \f$\theta/2\f$ and
         *        flag `fixed`. The default value of `fixed` is false.
         */
        QRotationGate2( const int* qubits , const angle_type& angle ,
                        const bool fixed = false )
        : QRotation< real_type >( angle , fixed )
        {
          setQubits( qubits ) ;
        } // QRotationGate2(qubits,angle,fixed)

        /**
         * \brief Constructs a 2-qubit rotation gate on the given qubits
         *        `qubits` with value `theta` = \f$\theta\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        QRotationGate2( const int* qubits , const real_type theta ,
                        const bool fixed = false )
        : QRotation< real_type >( theta , fixed )
        {
          setQubits( qubits ) ;
        } // QRotationGate2(qubits,theta,fixed)

        /**
         * \brief Constructs a 2-qubit rotation gate on the given qubits
         *        `qubits` with values `cos` = \f$\cos(\theta/2)\f$ and
         *        `sin` = \f$\sin(\theta/2)\f$, and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        QRotationGate2( const int* qubits , const real_type cos ,
                        const real_type sin , const bool fixed = false )
        : QRotation< real_type >( cos , sin , fixed )
        {
          setQubits( qubits ) ;
        } // QRotationGate2(qubits,cos,sin,fixed)

        /**
         * \brief Constructs a 2-qubit rotation gate on the given
         *        qubits `qubit0` and `qubit1` with quantum angle
         *        `angle` = \f$\theta/2\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        QRotationGate2( const int qubit0 , const int qubit1 ,
                        const angle_type& angle , const bool fixed = false )
        : QRotation< real_type >( angle , fixed )
        {
          const int qubits[] = { qubit0 , qubit1 } ;
          setQubits( &qubits[0] ) ;
        } // QRotationGate2(qubit0,qubit1,angle,fixed)

        /**
         * \brief Constructs a 2-qubit rotation gate on the given qubits
         *        `qubit0` and `qubit1` with value `theta` = \f$\theta\f$
         *        and flag `fixed`. The default value of `fixed` is false.
         */
        QRotationGate2( const int qubit0 , const int qubit1 ,
                        const real_type theta , const bool fixed = false )
        : QRotation< real_type >( theta , fixed )
        {
          const int qubits[] = { qubit0 , qubit1 } ;
          setQubits( &qubits[0] ) ;
        } // QRotationGate2(qubit0,qubit1,theta,fixed)

        /**
         * \brief Constructs a 2-qubit rotation gate on the given qubits
         *        `qubit0` and `qubit1` with values `cos` = \f$\cos(\theta/2)\f$
         *        and `sin` = \f$\sin(\theta/2)\f$, and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        QRotationGate2( const int qubit0 , const int qubit1 ,
                        const real_type cos , const real_type sin ,
                        const bool fixed = false )
        : QRotation< real_type >( cos , sin , fixed )
        {
          const int qubits[] = { qubit0 , qubit1 } ;
          setQubits( &qubits[0] ) ;
        } // QRotationGate2(qubit0,qubit1,cos,sin,fixed)

        // nbQubits

        // fixed
        inline bool fixed() const override {
          return QRotation< real_type >::fixed() ;
        }

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

        // apply

        // print

        // toQASM

        /// Checks if `other` is equal to this 2-qubit rotation gate.
        inline bool operator==( const QRotationGate2< T >& other ) const {
          return other.angle() == this->angle() ;
        }

        /// Checks if `other` is different from this 2-qubit rotation gate.
        inline bool operator!=( const QRotationGate2< T >& other ) const {
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

      protected:
        int  qubits_[2] ;  ///< Qubits of this 2-qubit rotation gate.

    } ; // class QRotationGate2

  } // namespace qgates

} // namespace qclab

#endif

