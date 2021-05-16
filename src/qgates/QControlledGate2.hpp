//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_QControlledGate2_hpp
#define qclab_qgates_QControlledGate2_hpp

#include "QGate2.hpp"
#include "QGate1.hpp"
#include "dense/kron.hpp"
#include <cassert>

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
        void apply( Side side , Op op , const int nbQubits ,
                    qclab::dense::SquareMatrix< T >& matrix ) const override {
          using matrix_type = qclab::dense::SquareMatrix< T > ;
          assert( nbQubits >= 2 ) ;
          assert( matrix.size() == 1 << nbQubits ) ;
          const auto qubits = this->qubits() ;
          assert( qubits[0] < nbQubits ) ; assert( qubits[1] < nbQubits ) ;
          if ( qubits[0] + 1 == qubits[1] ) {
            // nearest neighbor qubits
            QGate2< T >::apply( side , op , nbQubits , matrix ) ;
            return ;
          }
          // 2x2 matrices
          matrix_type  E0( 1 , 0 ,
                           0 , 0 ) ;
          matrix_type  E1( 0 , 0 ,
                           0 , 1 ) ;
          matrix_type  I2 = qclab::dense::eye< T >( 2 ) ;
          matrix_type  mat1 = this->gate()->matrix() ;
          // operation
          qclab::dense::operateInPlace( op , mat1 ) ;
          // linear combination of kronecker products
          const int s = qubits[1] - qubits[0] + 1 ;
          matrix_type  mats( 1 << s ) ;
          matrix_type  Imid = qclab::dense::eye< T >( 1 << ( s - 2 ) ) ;
          if ( controlState_ == 0 ) {
            if ( control() < target() ) {
              // kron( E0 , Imid , mat1 ) + kron( E1 , Imid , I2 )
              auto tmp0 = qclab::dense::kron( E0 , Imid ) ;
              auto tmp1 = qclab::dense::kron( E1 , Imid ) ;
              mats = qclab::dense::kron( tmp0 , mat1 ) +
                     qclab::dense::kron( tmp1 ,  I2  ) ;

            } else {
              // kron( mat1 , Imid , E0 ) + kron( I2 , Imid , E1 )
              auto tmp0 = qclab::dense::kron( mat1 , Imid ) ;
              auto tmp1 = qclab::dense::kron(  I2  , Imid ) ;
              mats = qclab::dense::kron( tmp0 , E0 ) +
                     qclab::dense::kron( tmp1 , E1 ) ;
            }
          } else {
            if ( control() < target() ) {
              // kron( E0 , Imid , I2 ) + kron( E1 , Imid , mat1 )
              auto tmp0 = qclab::dense::kron( E0 , Imid ) ;
              auto tmp1 = qclab::dense::kron( E1 , Imid ) ;
              mats = qclab::dense::kron( tmp0 ,  I2  ) +
                     qclab::dense::kron( tmp1 , mat1 ) ;
            } else {
              // kron( I2 , Imid , E0 ) + kron( mat1 , Imid , E1 )
              auto tmp0 = qclab::dense::kron(  I2  , Imid ) ;
              auto tmp1 = qclab::dense::kron( mat1 , Imid ) ;
              mats = qclab::dense::kron( tmp0 , E0 ) +
                     qclab::dense::kron( tmp1 , E1 ) ;
            }
          }
          // kron( Ileft , mats , Iright )
          qclab::dense::SquareMatrix< T >  matn( 1 << nbQubits ) ;
          if ( qubits[0] == 0 && qubits[1] == nbQubits - 1 ) {
            matn = mats ;
          } else if ( qubits[0] == 0 ) {
            auto Iright = qclab::dense::eye< T >( 1 << ( nbQubits - s ) ) ;
            qclab::dense::kron( mats , Iright , matn ) ;
          } else if ( qubits[1] == nbQubits - 1 ) {
            auto Ileft = qclab::dense::eye< T >( 1 << ( nbQubits - s ) ) ;
            qclab::dense::kron( Ileft , mats , matn ) ;
          } else {
            auto Ileft  = qclab::dense::eye< T >( 1 << qubits[0] ) ;
            auto Iright = qclab::dense::eye< T >( 1 << ( nbQubits-qubits[1]-1));
            auto tmp = qclab::dense::kron( Ileft , mats ) ;
            qclab::dense::kron( tmp , Iright , matn ) ;
          }
          // side
          if ( side == Side::Left ) {
            matrix *= matn ;
          } else {
            matrix = matn * matrix ;
          }
        }

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

#endif

