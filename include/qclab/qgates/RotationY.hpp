//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#pragma once

#include "qclab/qgates/QRotationGate1.hpp"

namespace qclab {

  namespace qgates {

    /**
     * \class RotationY
     * \brief 1-qubit rotation gate about the Y axis.
     */
    template <typename T>
    class RotationY : public QRotationGate1< T >
    {

      public:
        /// Real value type of this 1-qubit Y-rotation gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum rotation type of this 1-qubit Y-rotation gate.
        using rotation_type = qclab::QRotation< real_type > ;

        /**
         * \brief Default constructor. Constructs a 1-qubit Y-rotation gate on
         *        qubit 0 with parameter \f$\theta = 0\f$.
         */
        RotationY()
        : QRotationGate1< T >()
        { } // RotationY()

        /**
         * \brief Constructs a 1-qubit Y-rotation gate on qubit 0 with the given
         *        quantum rotation `rot` = \f$\theta\f$.
         */
        RotationY( const rotation_type& rot )
        : QRotationGate1< T >( 0 , rot )
        { } // RotationY(rot)

        /**
         * \brief Constructs a 1-qubit Y-rotation gate on qubit 0 with the given
         *        value `theta` = \f$\theta\f$.
         */
        RotationY( const real_type theta )
        : QRotationGate1< T >( 0 , theta )
        { } // RotationY(theta)

        /**
         * \brief Constructs a 1-qubit Y-rotation gate on qubit 0 with the given
         *        values `cos` = \f$\cos(\theta/2)\f$ and
         *               `sin` = \f$\sin(\theta/2)\f$.
         */
        RotationY( const real_type cos , const real_type sin )
        : QRotationGate1< T >( 0 , cos , sin )
        { } // RotationY(cos,sin)

        /**
         * \brief Constructs a 1-qubit Y-rotation gate on the given qubit
         *        `qubit` with quantum rotation `rot` = \f$\theta\f$ and
         *        flag `fixed`. The default value of `fixed` is false.
         */
        RotationY( const int qubit , const rotation_type& rot ,
                   const bool fixed = false )
        : QRotationGate1< T >( qubit , rot , fixed )
        { } // RotationY(qubit,rot,fixed)

        /**
         * \brief Constructs a 1-qubit Y-rotation gate on the given qubit
         *        `qubit` with value `theta` = \f$\theta\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationY( const int qubit , const real_type theta ,
                   const bool fixed = false )
        : QRotationGate1< T >( qubit , theta , fixed )
        { } // RotationY(qubit,theta,fixed)

        /**
         * \brief Constructs a 1-qubit Y-rotation gate on the given qubit
         *        `qubit` with values `cos` = \f$\cos(\theta/2)\f$ and
         *        `sin` = \f$\sin(\theta/2)\f$, and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationY( const int qubit , const real_type cos , const real_type sin ,
                   const bool fixed = false )
        : QRotationGate1< T >( qubit , cos , sin , fixed )
        { } // RotationY(qubit,cos,sin,fixed)

        // nbQubits

        // fixed

        // controlled

        // qubit

        // setQubit

        // qubits

        // setQubits

        // matrix
        qclab::dense::SquareMatrix< T > matrix() const override {
          return qclab::dense::SquareMatrix< T >( this->cos() , -this->sin() ,
                                                  this->sin() ,  this->cos() ) ;
        }

        // apply
        void apply( Op op , const int nbQubits , std::vector< T >& vector ,
                    const int offset = 0 ) const override ;

      #ifdef QCLAB_OMP_OFFLOADING
        // apply_device
        void apply_device( Op op , const int nbQubits , T* vector ,
                           const int offset = 0 ) const override ;
      #endif

        // apply
        void apply( Side side , Op op , const int nbQubits ,
                    qclab::dense::SquareMatrix< T >& matrix ,
                    const int offset = 0 ) const override ;

        // print

        // toQASM
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          stream << qasmRy( this->qubit_ + offset , this->theta() ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using Y = RotationY< T > ;
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

        /// Multiplies `rhs` to this 1-qubit Y-rotation gate.
        inline RotationY< T >& operator*=( const RotationY< T >& rhs ) {
          assert( this->qubit() == rhs.qubit() ) ;
          this->rotation_ *= rhs.rotation() ;
          return *this ;
        }

        /// Multiplies the inverse of `rhs` to this 1-qubit Y-rotation gate.
        inline RotationY< T >& operator/=( const RotationY< T >& rhs ) {
          assert( this->qubit() == rhs.qubit() ) ;
          this->rotation_ /= rhs.rotation() ;
          return *this ;
        }

        /// Multiplies `lhs` and `rhs`.
        friend RotationY< T > operator*( RotationY< T > lhs ,
                                         const RotationY< T >& rhs ) {
          assert( lhs.qubit() == rhs.qubit() ) ;
          lhs *= rhs ;
          return lhs ;
        }

        /// Multiplies `lhs` and the inverse of `rhs`.
        friend RotationY< T > operator/( RotationY< T > lhs ,
                                         const RotationY< T >& rhs ) {
          assert( lhs.qubit() == rhs.qubit() ) ;
          lhs /= rhs ;
          return lhs ;
        }

        /// Returns the inverse this 1-qubit Y-rotation gate.
        inline RotationY< T > inv() const {
          RotationY< T > gate( this->rotation().inv() ) ;
          return gate ;
        }

    } ; // class RotationY

  } // namespace qgates

} // namespace qclab

