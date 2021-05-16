//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_QGate2_hpp
#define qclab_qgates_QGate2_hpp

#include "QObject.hpp"
#include "dense/kron.hpp"
#include "dense/transpose.hpp"
#include <cassert>

namespace qclab {

  namespace qgates {

    /**
     * \class QGate2
     * \brief Base class for 2-qubit gates.
     */
    template <typename T>
    class QGate2 : public qclab::QObject< T >
    {

      public:
        // nbQubits
        inline int nbQubits() const override { return 2 ; }

        // fixed

        // controlled

        // qubit

        // setQubit
        inline void setQubit( const int qubit ) override { assert( false ) ; }

        // qubits

        // setQubits

        // matrix

        // apply
        void apply( Side side , Op op , const int nbQubits ,
                    qclab::dense::SquareMatrix< T >& matrix ) const override {
          assert( nbQubits >= 2 ) ;
          assert( matrix.size() == 1 << nbQubits ) ;
          const auto qubits = this->qubits() ;
          assert( qubits[0] < nbQubits ) ; assert( qubits[1] < nbQubits ) ;
          assert( qubits[0] + 1 == qubits[1] ) ;  // nearest neighbor qubtis
          // operation
          qclab::dense::SquareMatrix< T >  mat2 = this->matrix() ;
          qclab::dense::operateInPlace( op , mat2 ) ;
          // kron( Ileft , mat2 , Iright )
          qclab::dense::SquareMatrix< T >  matn( 1 << nbQubits ) ;
          if ( nbQubits == 2 ) {
            matn = mat2 ;
          } else if ( qubits[0] == 0 ) {
            auto Iright = qclab::dense::eye< T >( 1 << ( nbQubits - 2 ) ) ;
            qclab::dense::kron( mat2 , Iright , matn ) ;
          } else if ( qubits[1] == nbQubits - 1 ) {
            auto Ileft = qclab::dense::eye< T >( 1 << ( nbQubits - 2 ) ) ;
            qclab::dense::kron( Ileft , mat2 , matn ) ;
          } else {
            auto Ileft  = qclab::dense::eye< T >( 1 << qubits[0] ) ;
            auto Iright = qclab::dense::eye< T >( 1 << ( nbQubits-qubits[1]-1));
            auto tmp = qclab::dense::kron( Ileft , mat2 ) ;
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
          printMatrix4x4( this->matrix() ) ;
        }

        // toQASM

        // operator==

        // operator!=

        // equals

    } ; // class QGate2

  } // namespace qgates

} // namespace qclab

#endif

