//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_RotationXX_hpp
#define qclab_qgates_RotationXX_hpp

#include "QRotationGate2.hpp"
#include <cassert>

namespace qclab {

  namespace qgates {

    /**
     * \class RotationXX
     * \brief 2-qubit rotation gate about XX.
     */
    template <typename T>
    class RotationXX : public QRotationGate2< T >
    {

      public:
        /// Real value type of this 2-qubit X-rotation gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum angle type of this 2-qubit X-rotation gate.
        using angle_type = qclab::QAngle< real_type > ;

        /**
         * \brief Default constructor. Constructs a 2-qubit X-rotation gate on
         *        qubits 0 and 1 with parameter \f$\theta = 0\f$.
         */
        RotationXX()
        : QRotationGate2< T >()
        { } // RotationXX()

        /**
         * \brief Constructs a 2-qubit X-rotation gate on qubits 0 and 1 with
         *        the given quantum angle `angle` = \f$\theta/2\f$.
         */
        RotationXX( const angle_type& angle )
        : QRotationGate2< T >( 0 , 1 , angle )
        { } // RotationXX(angle)

        /**
         * \brief Constructs a 2-qubit X-rotation gate on qubits 0 and 1 with
         *        the given value `theta` = \f$\theta\f$.
         */
        RotationXX( const real_type theta )
        : QRotationGate2< T >( 0 , 1 , theta )
        { } // RotationXX(theta)

        /**
         * \brief Constructs a 2-qubit X-rotation gate on qubits 0 and 1 with
         *        the given values `cos` = \f$\cos(\theta/2)\f$ and
         *                         `sin` = \f$\sin(\theta/2)\f$.
         */
        RotationXX( const real_type cos , const real_type sin )
        : QRotationGate2< T >( 0 , 1 , cos , sin )
        { } // RotationXX(cos,sin)

        /**
         * \brief Constructs a 2-qubit X-rotation gate on the given qubits
         *        `qubits` with quantum angle `angle` = \f$\theta/2\f$ and
         *        flag `fixed`. The default value of `fixed` is false.
         */
        RotationXX( const int* qubits , const angle_type& angle ,
                    const bool fixed = false )
        : QRotationGate2< T >( qubits , angle , fixed )
        { } // RotationXX(qubits,angle,fixed)

        /**
         * \brief Constructs a 2-qubit X-rotation gate on the given qubits
         *        `qubits` with value `theta` = \f$\theta\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationXX( const int* qubits , const real_type theta ,
                    const bool fixed = false )
        : QRotationGate2< T >( qubits , theta , fixed )
        { } // RotationXX(qubits,theta,fixed)

        /**
         * \brief Constructs a 2-qubit X-rotation gate on the given qubits
         *        `qubits` with values `cos` = \f$\cos(\theta/2)\f$ and
         *        `sin` = \f$\sin(\theta/2)\f$, and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationXX( const int* qubits , const real_type cos ,
                    const real_type sin , const bool fixed = false )
        : QRotationGate2< T >( qubits , cos , sin , fixed )
        { } // RotationXX(qubits,cos,sin,fixed)

        /**
         * \brief Constructs a 2-qubit X-rotation gate on the given
         *        qubits `qubit0` and `qubit1` with quantum angle
         *        `angle` = \f$\theta/2\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationXX( const int qubit0 , const int qubit1 ,
                    const angle_type& angle , const bool fixed = false )
        : QRotationGate2< T >( qubit0 , qubit1 , angle , fixed )
        { } // RotationXX(qubit0,qubit1,angle,fixed)

        /**
         * \brief Constructs a 2-qubit X-rotation gate on the given qubits
         *        `qubit0` and `qubit1`  with value `theta` = \f$\theta\f$
         *        and flag `fixed`. The default value of `fixed` is false.
         */
        RotationXX( const int qubit0 , const int qubit1 ,
                    const real_type theta , const bool fixed = false )
        : QRotationGate2< T >( qubit0 , qubit1 , theta , fixed )
        { } // RotationXX(qubit0,qubit1,theta,fixed)

        /**
         * \brief Constructs a 2-qubit X-rotation gate on the given qubits
         *        `qubit0` and `qubit1` with values `cos` = \f$\cos(\theta/2)\f$
         *        and `sin` = \f$\sin(\theta/2)\f$, and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationXX( const int qubit0 , const int qubit1 , const real_type cos ,
                    const real_type sin , const bool fixed = false )
        : QRotationGate2< T >( qubit0 , qubit1 , cos , sin , fixed )
        { } // RotationXX(qubit0,qubit1,cos,sin,fixed)

        // nbQubits

        // fixed

        // controlled

        // qubit

        // setQubit

        // qubits

        // setQubits

        // matrix
        qclab::dense::SquareMatrix< T > matrix() const override {
          const T d = this->cos() ;
          const T o = T(0,-this->sin()) ;
          return qclab::dense::SquareMatrix< T >( d , 0 , 0 , o ,
                                                  0 , d , o , 0 ,
                                                  0 , o , d , 0 ,
                                                  o , 0 , 0 , d ) ;
        }

        // apply

        // print

        // toQASM
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          auto qubits = this->qubits() ;
          stream << qasmRxx( qubits[0] + offset ,
                             qubits[1] + offset , this->theta() ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using X = RotationXX< T > ;
          if ( const X* p = dynamic_cast< const X* >( &other ) ) {
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

        /// Multiplies `rhs` to this 2-qubit X-rotation gate.
        inline RotationXX< T >& operator*=( const RotationXX< T >& rhs ) {
          this->update( this->angle() + rhs.angle() ) ;
          return *this ;
        }

        /// Multiplies the inverse of `rhs` to this 2-qubit X-rotation gate.
        inline RotationXX< T >& operator/=( const RotationXX< T >& rhs ) {
          this->update( this->angle() - rhs.angle() ) ;
          return *this ;
        }

        /// Multiplies `lhs` and `rhs`.
        friend RotationXX< T > operator*( RotationXX< T > lhs ,
                                         const RotationXX< T >& rhs ) {
          lhs *= rhs ;
          return lhs ;
        }

        /// Multiplies `lhs` and the inverse of `rhs`.
        friend RotationXX< T > operator/( RotationXX< T > lhs ,
                                         const RotationXX< T >& rhs ) {
          lhs /= rhs ;
          return lhs ;
        }

        /// Returns the inverse this 2-qubit X-rotation gate.
        inline RotationXX< T > inv() const {
          RotationXX< T > rotation( -this->angle() ) ;
          return rotation ;
        }

    } ; // class RotationXX

  } // namespace qgates

} // namespace qclab

#endif

