//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#pragma once

#include "qclab/qgates/QGate1.hpp"
#include "qclab/qgates/QGate2.hpp"
#include "qclab/dense/kron.hpp"

namespace qclab {

  namespace qgates {

    /**
     * \class QControlledGate2
     * \brief Base class for 2-qubit gates of controlled 1-qubit gates.
     */
    template <typename T>
    class QControlledGate2 : public QGate2< T >
    {

      public:
        /// Gate type of the controlled 1-qubit gate of this 2-qubit gate.
        using gate_type = QGate1< T > ;

        /**
         * \brief Constructs a 2-qubit gate with the given control qubit `qubit`
         *        on the control state `controlState`.
         *        The default control state is 1.
         */
        QControlledGate2( const int control , const int controlState = 1 )
        : control_( control )
        , controlState_( controlState )
        {
          assert( control >= 0 ) ;
          assert( ( controlState == 0 ) || ( controlState == 1 ) ) ;
        } // QControlledGate2(control,controlState)

        // nbQubits

        // fixed

        // controlled
        inline bool controlled() const override { return true ; }

        // qubit

        // setQubit

        // qubits
        std::vector< int > qubits() const override {
          return std::vector< int >( { std::min( control() , this->target() ) ,
                                       std::max( control() , this->target() )});
        }

        // setQubits
        inline void setQubits( const int* qubits ) override {
          assert( qubits[0] >= 0 ) ; assert( qubits[1] >= 0 ) ;
          assert( qubits[0] != qubits[1] ) ;
          control_ = qubits[0] ;
          this->setTarget( qubits[1] ) ;
        }

        // matrix
        qclab::dense::SquareMatrix< T > matrix() const override {
          using matrix_type = qclab::dense::SquareMatrix< T > ;
          const matrix_type  I2 = qclab::dense::eye< T >( 2 ) ;
          const matrix_type  E0( 1 , 0 ,
                                 0 , 0 ) ;
          const matrix_type  E1( 0 , 0 ,
                                 0 , 1 ) ;
          const matrix_type  CG = this->gate()->matrix() ;
          if ( controlState_ == 0 ) {
            if ( control() < target() ) {
              return qclab::dense::kron( E0 , CG ) +
                     qclab::dense::kron( E1 , I2 ) ;
            } else {
              return qclab::dense::kron( CG , E0 ) +
                     qclab::dense::kron( I2 , E1 ) ;
            }
          } else {
            if ( control() < target() ) {
              return qclab::dense::kron( E0 , I2 ) +
                     qclab::dense::kron( E1 , CG ) ;
            } else {
              return qclab::dense::kron( I2 , E0 ) +
                     qclab::dense::kron( CG , E1 ) ;
            }
          }
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

        // operator==

        // operator!=

        // equals

        /// Returns the control qubit of this 2-qubit gate.
        inline int control() const { return control_ ; }

        /// Returns the target qubit of this 2-qubit gate.
        virtual int target() const = 0 ;

        /// Returns the control state of this 2-qubit gate.
        inline int controlState() const { return controlState_ ; }

        /**
         * \brief Returns a const pointer to the controlled 1-qubit gate of this
         *        2-qubit gate.
         */
        virtual const gate_type* gate() const = 0 ;

        /// Sets the control of this controlled gate to the given `control`.
        void setControl( const int control ) {
          assert( control >= 0 ) ; assert( control != target() ) ;
          control_ = control ;
        }

        /// Sets the target of this controlled gate to the given `target`.
        virtual void setTarget( const int target ) = 0 ;

        /// Sets the control state of this controlled gate to `controlState`.
        void setControlState( const int controlState ) {
          assert( ( controlState == 0 ) || ( controlState == 1 ) ) ;
          controlState_ = controlState ;
        }

      protected:
        int  control_ ;       ///< Control qubit of this 2-qubit gate.
        int  controlState_ ;  ///< Control state of this 2-qubit gate.

    } ; // class QControlledGate2

  } // namespace qgates

} // namespace qclab

