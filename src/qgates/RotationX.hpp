//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_RotationX_hpp
#define qclab_qgates_RotationX_hpp

#include "QRotationGate1.hpp"
#include <cassert>

namespace qclab {

  namespace qgates {

    /**
     * \class RotationX
     * \brief 1-qubit rotation gate about the X axis.
     */
    template <typename T>
    class RotationX : public QRotationGate1< T >
    {

      public:
        /// Real value type of this 1-qubit X-rotation gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum angle type of this 1-qubit X-rotation gate.
        using angle_type = qclab::QAngle< real_type > ;

        /**
         * \brief Default constructor. Constructs a 1-qubit X-rotation gate on
         *        qubit 0 with parameter \f$\theta = 0\f$.
         */
        RotationX()
        : QRotationGate1< T >()
        { } // RotationX()

        /**
         * \brief Constructs a 1-qubit X-rotation gate on qubit 0 with the given
         *        quantum angle `angle` = \f$\theta/2\f$.
         */
        RotationX( const angle_type& angle )
        : QRotationGate1< T >( 0 , angle )
        { } // RotationX(angle)

        /**
         * \brief Constructs a 1-qubit X-rotation gate on qubit 0 with the given
         *        value `theta` = \f$\theta\f$.
         */
        RotationX( const real_type theta )
        : QRotationGate1< T >( 0 , theta )
        { } // RotationX(theta)

        /**
         * \brief Constructs a 1-qubit X-rotation gate on qubit 0 with the given
         *        values `cos` = \f$\cos(\theta/2)\f$ and
         *               `sin` = \f$\sin(\theta/2)\f$.
         */
        RotationX( const real_type cos , const real_type sin )
        : QRotationGate1< T >( 0 , cos , sin )
        { } // RotationX(cos,sin)

        /**
         * \brief Constructs a 1-qubit X-rotation gate on the given qubit
         *        `qubit` with quantum angle `angle` = \f$\theta/2\f$ and
         *        flag `fixed`. The default value of `fixed` is false.
         */
        RotationX( const int qubit , const angle_type& angle ,
                   const bool fixed = false )
        : QRotationGate1< T >( qubit , angle , fixed )
        { } // RotationX(qubit,angle,fixed)

        /**
         * \brief Constructs a 1-qubit X-rotation gate on the given qubit
         *        `qubit` with value `theta` = \f$\theta\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationX( const int qubit , const real_type theta ,
                   const bool fixed = false )
        : QRotationGate1< T >( qubit , theta , fixed )
        { } // RotationX(qubit,theta,fixed)

        /**
         * \brief Constructs a 1-qubit X-rotation gate on the given qubit
         *        `qubit` with values `cos` = \f$\cos(\theta/2)\f$ and
         *        `sin` = \f$\sin(\theta/2)\f$, and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationX( const int qubit , const real_type cos , const real_type sin ,
                   const bool fixed = false )
        : QRotationGate1< T >( qubit , cos , sin , fixed )
        { } // RotationX(qubit,cos,sin,fixed)

        // nbQubits

        // fixed

        // controlled

        // qubit

        // setQubit

        // qubits

        // setQubits

        // matrix
        qclab::dense::SquareMatrix< T > matrix() const override {
          using SqMat = qclab::dense::SquareMatrix< T > ;
          return SqMat(      this->cos()  , T(0,-this->sin()) ,
                        T(0,-this->sin()) ,      this->cos()  ) ;
        }

        // apply

        // print

        // toQASM
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          stream << qasmRx( this->qubit_ + offset , this->theta() ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using X = RotationX< T > ;
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

    } ; // class RotationX

  } // namespace qgates

} // namespace qclab

#endif

