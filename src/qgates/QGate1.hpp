//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_QGate1_hpp
#define qclab_qgates_QGate1_hpp

#include "QObject.hpp"
#include "dense/kron.hpp"
#include "dense/transpose.hpp"
#include <cassert>

namespace qclab {

  /**
   * Namespace qclab::qgates.
   */
  namespace qgates {

    /**
     * \class QGate1
     * \brief Base class for 1-qubit gates.
     */
    template <typename T>
    class QGate1 : public qclab::QObject< T >
    {

      public:
        /// Default contructor. Constructs a 1-qubit gate on qubit 0.
        QGate1()
        : qubit_( 0 )
        { } // QGate1()

        /// Constructs a 1-qubit gate on the given qubit `qubit`.
        QGate1( const int qubit )
        : qubit_( qubit )
        {
          assert( qubit >= 0 ) ;
        } // QGate1(qubit)

        // nbQubits
        inline int nbQubits() const override { return 1 ; }

        // fixed

        // controlled
        inline bool controlled() const override { return false ; }

        // qubit
        inline int qubit() const override { return qubit_ ; }

        // setQubit
        inline void setQubit( const int qubit ) override {
          assert( qubit >= 0 ) ;
          qubit_ = qubit ;
        }

        // qubits
        std::vector< int > qubits() const override {
          return std::vector< int >( { qubit_ } ) ;
        }

        // setQubits
        inline void setQubits( const int* qubits ) override {
          assert( qubits[0] >= 0 ) ;
          qubit_ = qubits[0] ;
        }

        // matrix

        // apply
        void apply( Side side , Op op , const int nbQubits ,
                    qclab::dense::SquareMatrix< T >& matrix ) const override {
          assert( nbQubits >= 1 ) ;
          assert( matrix.size() == 1 << nbQubits ) ;
          const int qubit = this->qubit() ;  assert( qubit < nbQubits ) ;
          // operation
          qclab::dense::SquareMatrix< T >  mat1 = this->matrix() ;
          qclab::dense::operateInPlace( op , mat1 ) ;
          // kron( Ileft , mat1 , Iright )
          qclab::dense::SquareMatrix< T >  matn( 1 << nbQubits ) ;
          if ( nbQubits == 1 ) {
            matn = mat1 ;
          } else if ( qubit == 0 ) {
            auto Iright = qclab::dense::eye< T >( 1 << ( nbQubits - 1 ) ) ;
            qclab::dense::kron( mat1 , Iright , matn ) ;
          } else if ( qubit == nbQubits - 1 ) {
            auto Ileft = qclab::dense::eye< T >( 1 << ( nbQubits - 1 ) ) ;
            qclab::dense::kron( Ileft , mat1 , matn ) ;
          } else {
            auto Ileft  = qclab::dense::eye< T >( 1 << qubit ) ;
            auto Iright = qclab::dense::eye< T >( 1 << ( nbQubits - qubit - 1));
            auto tmp = qclab::dense::kron( Ileft , mat1 ) ;
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
        void print() const override {
          printMatrix2x2( this->matrix() ) ;
        }

        // toQASM

        // operator==

        // operator!=

        // equals

      protected:
        int  qubit_ ;  ///< Qubit of this 1-qubit gate.

    } ; // class QGate1

  } // namespace qgates

} // namespace qclab

#endif

