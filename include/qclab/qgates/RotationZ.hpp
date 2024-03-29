//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#pragma once

#include "qclab/qgates/QRotationGate1.hpp"

namespace qclab {

  namespace qgates {

    /**
     * \class RotationZ
     * \brief 1-qubit rotation gate about the Z axis.
     */
    template <typename T>
    class RotationZ : public QRotationGate1< T >
    {

      public:
        /// Real value type of this 1-qubit Z-rotation gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum rotation type of this 1-qubit Z-rotation gate.
        using rotation_type = qclab::QRotation< real_type > ;

        /**
         * \brief Default constructor. Constructs a 1-qubit Z-rotation gate on
         *        qubit 0 with parameter \f$\theta = 0\f$.
         */
        RotationZ()
        : QRotationGate1< T >()
        { } // RotationZ()

        /**
         * \brief Constructs a 1-qubit Z-rotation gate on qubit 0 with the given
         *        quantum rotation `rot` = \f$\theta\f$.
         */
        RotationZ( const rotation_type& rot )
        : QRotationGate1< T >( 0 , rot )
        { } // RotationZ(rot)

        /**
         * \brief Constructs a 1-qubit Z-rotation gate on qubit 0 with the given
         *        value `theta` = \f$\theta\f$.
         */
        RotationZ( const real_type theta )
        : QRotationGate1< T >( 0 , theta )
        { } // RotationZ(theta)

        /**
         * \brief Constructs a 1-qubit Z-rotation gate on qubit 0 with the given
         *        values `cos` = \f$\cos(\theta/2)\f$ and
         *               `sin` = \f$\sin(\theta/2)\f$.
         */
        RotationZ( const real_type cos , const real_type sin )
        : QRotationGate1< T >( 0 , cos , sin )
        { } // RotationZ(cos,sin)

        /**
         * \brief Constructs a 1-qubit Z-rotation gate on the given qubit
         *        `qubit` with quantum rotation `rot` = \f$\theta\f$ and
         *        flag `fixed`. The default value of `fixed` is false.
         */
        RotationZ( const int qubit , const rotation_type& rot ,
                   const bool fixed = false )
        : QRotationGate1< T >( qubit , rot , fixed )
        { } // RotationZ(qubit,rot,fixed)

        /**
         * \brief Constructs a 1-qubit Z-rotation gate on the given qubit
         *        `qubit` with value `theta` = \f$\theta\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationZ( const int qubit , const real_type theta ,
                   const bool fixed = false )
        : QRotationGate1< T >( qubit , theta , fixed )
        { } // RotationZ(qubit,theta,fixed)

        /**
         * \brief Constructs a 1-qubit Z-rotation gate on the given qubit
         *        `qubit` with values `cos` = \f$\cos(\theta/2)\f$ and
         *        `sin` = \f$\sin(\theta/2)\f$, and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        RotationZ( const int qubit , const real_type cos , const real_type sin ,
                   const bool fixed = false )
        : QRotationGate1< T >( qubit , cos , sin , fixed )
        { } // RotationZ(qubit,cos,sin,fixed)

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
          return SqMat( T(this->cos(),-this->sin()) , T(0) ,
                        T(0) , T(this->cos(), this->sin()) ) ;
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
          stream << qasmRz( this->qubit_ + offset , this->theta() ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using Z = RotationZ< T > ;
          if ( const Z* p = dynamic_cast< const Z* >( &other ) ) {
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

        /// Multiplies `rhs` to this 1-qubit Z-rotation gate.
        inline RotationZ< T >& operator*=( const RotationZ< T >& rhs ) {
          assert( this->qubit() == rhs.qubit() ) ;
          this->rotation_ *= rhs.rotation() ;
          return *this ;
        }

        /// Multiplies the inverse of `rhs` to this 1-qubit Z-rotation gate.
        inline RotationZ< T >& operator/=( const RotationZ< T >& rhs ) {
          assert( this->qubit() == rhs.qubit() ) ;
          this->rotation_ /= rhs.rotation() ;
          return *this ;
        }

        /// Multiplies `lhs` and `rhs`.
        friend RotationZ< T > operator*( RotationZ< T > lhs ,
                                         const RotationZ< T >& rhs ) {
          assert( lhs.qubit() == rhs.qubit() ) ;
          lhs *= rhs ;
          return lhs ;
        }

        /// Multiplies `lhs` and the inverse of `rhs`.
        friend RotationZ< T > operator/( RotationZ< T > lhs ,
                                         const RotationZ< T >& rhs ) {
          assert( lhs.qubit() == rhs.qubit() ) ;
          lhs /= rhs ;
          return lhs ;
        }

        /// Returns the inverse this 1-qubit Z-rotation gate.
        inline RotationZ< T > inv() const {
          RotationZ< T > gate( this->rotation().inv() ) ;
          return gate ;
        }

    } ; // class RotationZ

  } // namespace qgates

} // namespace qclab

