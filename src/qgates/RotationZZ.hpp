//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_RotationZZ_hpp
#define qclab_qgates_RotationZZ_hpp

#include "QRotationGate2.hpp"
#include <cassert>

namespace qclab {

  namespace qgates {

    /**
     * \class RotationZZ
     * \brief 2-qubit rotation gate about ZZ.
     */
    template <typename T>
    class RotationZZ : public QRotationGate2< T >
    {

      public:
        /// Real value type of this 2-qubit Z-rotation gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum angle type of this 2-qubit Z-rotation gate.
        using angle_type = qclab::QAngle< real_type > ;

        /**
         * \brief Default constructor. Constructs a 2-qubit Z-rotation gate on
         *        qubits 0 and 1 with parameter \f$\theta = 0\f$.
         */
        RotationZZ()
        : QRotationGate2< T >()
        { } // RotationZZ()

        /**
         * \brief Constructs a 2-qubit Z-rotation gate on qubits 0 and 1 with
         *        the given quantum angle `angle` = \f$\theta/2\f$.
         */
        RotationZZ( const angle_type& angle )
        : QRotationGate2< T >( 0 , 1 , angle )
        { } // RotationZZ(angle)

        /**
         * \brief Constructs a 2-qubit Z-rotation gate on qubits 0 and 1 with
         *        the given value `theta` = \f$\theta\f$.
         */
        RotationZZ( const real_type theta )
        : QRotationGate2< T >( 0 , 1 , theta )
        { } // RotationZZ(theta)

        /**
         * \brief Constructs a 2-qubit Z-rotation gate on qubits 0 and 1 with
         *        the given values `cos` = \f$\cos(\theta/2)\f$ and
         *                         `sin` = \f$\sin(\theta/2)\f$.
         */
        RotationZZ( const real_type cos , const real_type sin )
        : QRotationGate2< T >( 0 , 1 , cos , sin )
        { } // RotationZZ(cos,sin)

        /**
         * \brief Constructs a 2-qubit Z-rotation gate on the given qubits
         *        `qubits` with quantum angle `angle` = \f$\theta/2\f$ and
         *        flag `fixed`. The default value of `fixed` is false.
         */
        RotationZZ( const int* qubits , const angle_type& angle ,
                    const bool fixed = false )
        : QRotationGate2< T >( qubits , angle , fixed )
        { } // RotationZZ(qubits,angle,fixed)

        /**
         * \brief Constructs a 2-qubit Z-rotation gate on the given qubits
         *        `qubits` with value `theta` = \f$\theta\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationZZ( const int* qubits , const real_type theta ,
                    const bool fixed = false )
        : QRotationGate2< T >( qubits , theta , fixed )
        { } // RotationZZ(qubits,theta,fixed)

        /**
         * \brief Constructs a 2-qubit Z-rotation gate on the given qubits
         *        `qubits` with values `cos` = \f$\cos(\theta/2)\f$ and
         *        `sin` = \f$\sin(\theta/2)\f$, and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationZZ( const int* qubits , const real_type cos ,
                    const real_type sin , const bool fixed = false )
        : QRotationGate2< T >( qubits , cos , sin , fixed )
        { } // RotationZZ(qubits,cos,sin,fixed)

        /**
         * \brief Constructs a 2-qubit Z-rotation gate on the given
         *        qubits `qubit0` and `qubit1` with quantum angle
         *        `angle` = \f$\theta/2\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationZZ( const int qubit0 , const int qubit1 ,
                    const angle_type& angle , const bool fixed = false )
        : QRotationGate2< T >( qubit0 , qubit1 , angle , fixed )
        { } // RotationZZ(qubit0,qubit1,angle,fixed)

        /**
         * \brief Constructs a 2-qubit Z-rotation gate on the given qubits
         *        `qubit0` and `qubit1`  with value `theta` = \f$\theta\f$
         *        and flag `fixed`. The default value of `fixed` is false.
         */
        RotationZZ( const int qubit0 , const int qubit1 ,
                    const real_type theta , const bool fixed = false )
        : QRotationGate2< T >( qubit0 , qubit1 , theta , fixed )
        { } // RotationZZ(qubit0,qubit1,theta,fixed)

        /**
         * \brief Constructs a 2-qubit Z-rotation gate on the given qubits
         *        `qubit0` and `qubit1` with values `cos` = \f$\cos(\theta/2)\f$
         *        and `sin` = \f$\sin(\theta/2)\f$, and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationZZ( const int qubit0 , const int qubit1 , const real_type cos ,
                    const real_type sin , const bool fixed = false )
        : QRotationGate2< T >( qubit0 , qubit1 , cos , sin , fixed )
        { } // RotationZZ(qubit0,qubit1,cos,sin,fixed)

        // nbQubits

        // fixed

        // controlled

        // qubit

        // setQubit

        // qubits

        // setQubits

        // matrix
        qclab::dense::SquareMatrix< T > matrix() const override {
          const T a = T(this->cos(),-this->sin()) ;
          const T b = T(this->cos(), this->sin()) ;
          return qclab::dense::SquareMatrix< T >( a , 0 , 0 , 0 ,
                                                  0 , b , 0 , 0 ,
                                                  0 , 0 , b , 0 ,
                                                  0 , 0 , 0 , a ) ;
        }

        // apply

        // print

        // toQASM
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          auto qubits = this->qubits() ;
          stream << qasmRzz( qubits[0] + offset ,
                             qubits[1] + offset , this->theta() ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using Z = RotationZZ< T > ;
          if ( const Z* p = dynamic_cast< const Z* >( &other ) ) {
            return ( p->angle() == this->angle() ) ;
          }
          return false ;
        }

        // angle

        // theta

        // cos

        // sin

        // update(angle)

        // update(theta)

        // update(cos,sin)

    } ; // class RotationZZ

  } // namespace qgates

} // namespace qclab

#endif

