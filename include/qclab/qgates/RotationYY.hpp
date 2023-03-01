//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#pragma once

#include "qclab/qgates/QRotationGate2.hpp"

namespace qclab {

  namespace qgates {

    /**
     * \class RotationYY
     * \brief 2-qubit rotation gate about YY.
     */
    template <typename T>
    class RotationYY : public QRotationGate2< T >
    {

      public:
        /// Real value type of this 2-qubit Y-rotation gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum rotation type of this 2-qubit Y-rotation gate.
        using rotation_type = qclab::QRotation< real_type > ;

        /**
         * \brief Default constructor. Constructs a 2-qubit Y-rotation gate on
         *        qubits 0 and 1 with parameter \f$\theta = 0\f$.
         */
        RotationYY()
        : QRotationGate2< T >()
        { } // RotationYY()

        /**
         * \brief Constructs a 2-qubit Y-rotation gate on qubits 0 and 1 with
         *        the given quantum rotation `rot` = \f$\theta\f$.
         */
        RotationYY( const rotation_type& rot )
        : QRotationGate2< T >( 0 , 1 , rot )
        { } // RotationYY(rot)

        /**
         * \brief Constructs a 2-qubit Y-rotation gate on qubits 0 and 1 with
         *        the given value `theta` = \f$\theta\f$.
         */
        RotationYY( const real_type theta )
        : QRotationGate2< T >( 0 , 1 , theta )
        { } // RotationYY(theta)

        /**
         * \brief Constructs a 2-qubit Y-rotation gate on qubits 0 and 1 with
         *        the given values `cos` = \f$\cos(\theta/2)\f$ and
         *                         `sin` = \f$\sin(\theta/2)\f$.
         */
        RotationYY( const real_type cos , const real_type sin )
        : QRotationGate2< T >( 0 , 1 , cos , sin )
        { } // RotationYY(cos,sin)

        /**
         * \brief Constructs a 2-qubit Y-rotation gate on the given qubits
         *        `qubits` with quantum rotation `rot` = \f$\theta\f$ and
         *        flag `fixed`. The default value of `fixed` is false.
         */
        RotationYY( const int* qubits , const rotation_type& rot ,
                    const bool fixed = false )
        : QRotationGate2< T >( qubits , rot , fixed )
        { } // RotationYY(qubits,rot,fixed)

        /**
         * \brief Constructs a 2-qubit Y-rotation gate on the given qubits
         *        `qubits` with value `theta` = \f$\theta\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationYY( const int* qubits , const real_type theta ,
                    const bool fixed = false )
        : QRotationGate2< T >( qubits , theta , fixed )
        { } // RotationYY(qubits,theta,fixed)

        /**
         * \brief Constructs a 2-qubit Y-rotation gate on the given qubits
         *        `qubits` with values `cos` = \f$\cos(\theta/2)\f$ and
         *        `sin` = \f$\sin(\theta/2)\f$, and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationYY( const int* qubits , const real_type cos ,
                    const real_type sin , const bool fixed = false )
        : QRotationGate2< T >( qubits , cos , sin , fixed )
        { } // RotationYY(qubits,cos,sin,fixed)

        /**
         * \brief Constructs a 2-qubit Y-rotation gate on the given
         *        qubits `qubit0` and `qubit1` with quantum rotation
         *        `rot` = \f$\theta\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationYY( const int qubit0 , const int qubit1 ,
                    const rotation_type& rot , const bool fixed = false )
        : QRotationGate2< T >( qubit0 , qubit1 , rot , fixed )
        { } // RotationYY(qubit0,qubit1,rot,fixed)

        /**
         * \brief Constructs a 2-qubit Y-rotation gate on the given qubits
         *        `qubit0` and `qubit1`  with value `theta` = \f$\theta\f$
         *        and flag `fixed`. The default value of `fixed` is false.
         */
        RotationYY( const int qubit0 , const int qubit1 ,
                    const real_type theta , const bool fixed = false )
        : QRotationGate2< T >( qubit0 , qubit1 , theta , fixed )
        { } // RotationYY(qubit0,qubit1,theta,fixed)

        /**
         * \brief Constructs a 2-qubit Y-rotation gate on the given qubits
         *        `qubit0` and `qubit1` with values `cos` = \f$\cos(\theta/2)\f$
         *        and `sin` = \f$\sin(\theta/2)\f$, and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationYY( const int qubit0 , const int qubit1 , const real_type cos ,
                    const real_type sin , const bool fixed = false )
        : QRotationGate2< T >( qubit0 , qubit1 , cos , sin , fixed )
        { } // RotationYY(qubit0,qubit1,cos,sin,fixed)

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
          const T o = T(0,this->sin()) ;
          return qclab::dense::SquareMatrix< T >( d , 0 , 0 , o ,
                                                  0 , d ,-o , 0 ,
                                                  0 ,-o , d , 0 ,
                                                  o , 0 , 0 , d ) ;
        }

        // apply

        // print

        // toQASM
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          auto qubits = this->qubits() ;
          stream << qasmRyy( qubits[0] + offset ,
                             qubits[1] + offset , this->theta() ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using Y = RotationYY< T > ;
          if ( const Y* p = dynamic_cast< const Y* >( &other ) ) {
            return ( p->rotation() == this->rotation() ) ;
          }
          return false ;
        }

        // rotation

        // theta

        // cos

        // sin

        // update(rot)

        // update(theta)

        // update(cos,sin)

        /// Multiplies `rhs` to this 2-qubit Y-rotation gate.
        inline RotationYY< T >& operator*=( const RotationYY< T >& rhs ) {
          assert( this->qubit() == rhs.qubit() ) ;
          this->rotation_ *= rhs.rotation() ;
          return *this ;
        }

        /// Multiplies the inverse of `rhs` to this 2-qubit Y-rotation gate.
        inline RotationYY< T >& operator/=( const RotationYY< T >& rhs ) {
          assert( this->qubit() == rhs.qubit() ) ;
          this->rotation_ /= rhs.rotation() ;
          return *this ;
        }

        /// Multiplies `lhs` and `rhs`.
        friend RotationYY< T > operator*( RotationYY< T > lhs ,
                                         const RotationYY< T >& rhs ) {
          assert( lhs.qubit() == rhs.qubit() ) ;
          lhs *= rhs ;
          return lhs ;
        }

        /// Multiplies `lhs` and the inverse of `rhs`.
        friend RotationYY< T > operator/( RotationYY< T > lhs ,
                                         const RotationYY< T >& rhs ) {
          assert( lhs.qubit() == rhs.qubit() ) ;
          lhs /= rhs ;
          return lhs ;
        }

        /// Returns the inverse this 2-qubit Y-rotation gate.
        inline RotationYY< T > inv() const {
          RotationYY< T > gate( this->rotation().inv() ) ;
          return gate ;
        }

    } ; // class RotationYY

  } // namespace qgates

} // namespace qclab

