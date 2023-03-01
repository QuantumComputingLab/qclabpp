//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#pragma once

#include "qclab/qgates/QRotationGate1.hpp"

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
        /// Quantum rotation type of this 1-qubit X-rotation gate.
        using rotation_type = qclab::QRotation< real_type > ;

        /**
         * \brief Default constructor. Constructs a 1-qubit X-rotation gate on
         *        qubit 0 with parameter \f$\theta = 0\f$.
         */
        RotationX()
        : QRotationGate1< T >()
        { } // RotationX()

        /**
         * \brief Constructs a 1-qubit X-rotation gate on qubit 0 with the given
         *        quantum rotation `rot` = \f$\theta\f$.
         */
        RotationX( const rotation_type& rot )
        : QRotationGate1< T >( 0 , rot )
        { } // RotationX(rot)

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
         *        `qubit` with quantum rotation `rot` = \f$\theta\f$ and
         *        flag `fixed`. The default value of `fixed` is false.
         */
        RotationX( const int qubit , const rotation_type& rot ,
                   const bool fixed = false )
        : QRotationGate1< T >( qubit , rot , fixed )
        { } // RotationX(qubit,rot,fixed)

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
          stream << qasmRx( this->qubit_ + offset , this->theta() ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using X = RotationX< T > ;
          if ( const X* p = dynamic_cast< const X* >( &other ) ) {
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

        /// Multiplies `rhs` to this 1-qubit X-rotation gate.
        inline RotationX< T >& operator*=( const RotationX< T >& rhs ) {
          assert( this->qubit() == rhs.qubit() ) ;
          this->rotation_ *= rhs.rotation() ;
          return *this ;
        }

        /// Multiplies the inverse of `rhs` to this 1-qubit X-rotation gate.
        inline RotationX< T >& operator/=( const RotationX< T >& rhs ) {
          assert( this->qubit() == rhs.qubit() ) ;
          this->rotation_ /= rhs.rotation() ;
          return *this ;
        }

        /// Multiplies `lhs` and `rhs`.
        friend RotationX< T > operator*( RotationX< T > lhs ,
                                         const RotationX< T >& rhs ) {
          assert( lhs.qubit() == rhs.qubit() ) ;
          lhs *= rhs ;
          return lhs ;
        }

        /// Multiplies `lhs` and the inverse of `rhs`.
        friend RotationX< T > operator/( RotationX< T > lhs ,
                                         const RotationX< T >& rhs ) {
          assert( lhs.qubit() == rhs.qubit() ) ;
          lhs /= rhs ;
          return lhs ;
        }

        /// Returns the inverse this 1-qubit X-rotation gate.
        inline RotationX< T > inv() const {
          RotationX< T > gate( this->rotation().inv() ) ;
          return gate ;
        }

    } ; // class RotationX

  } // namespace qgates

} // namespace qclab

