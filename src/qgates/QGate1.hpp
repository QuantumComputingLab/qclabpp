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
                    qclab::dense::SquareMatrix< T >& matrix ,
                    const int offset = 0 ) const override {
          assert( nbQubits >= 1 ) ;
          assert( matrix.size() == 1 << nbQubits ) ;
          const int qubit = this->qubit() + offset ;
          assert( qubit < nbQubits ) ;
          // operation
          qclab::dense::SquareMatrix< T >  mat1 = this->matrix() ;
          qclab::dense::operateInPlace( op , mat1 ) ;
          // side
          if ( side == Side::Left ) {
            if ( nbQubits == 1 ) {
              matrix *= mat1 ;
            } else {
              // matrix *= kron( Ileft , mat1 , Iright )
              const int64_t nLeft  = 1 << qubit ;
              const int64_t nRight = 1 << ( nbQubits - qubit - 1 ) ;
              #pragma omp parallel for
              for ( int64_t i = 0; i < matrix.rows(); i++ ) {
                for ( int64_t k = 0; k < nLeft; k++ ) {
                  for ( int64_t j = 0; j < nRight; j++ ) {
                    const int64_t j1 = j + 2 * k * nRight ;
                    const int64_t j2 = j1 + nRight ;
                    const T x1 = matrix(i,j1) ;
                    const T x2 = matrix(i,j2) ;
                    matrix(i,j1) = mat1(0,0) * x1 + mat1(1,0) * x2 ;
                    matrix(i,j2) = mat1(0,1) * x1 + mat1(1,1) * x2 ;
                  }
                }
              }
            }
          } else {
            if ( nbQubits == 1 ) {
              matrix = mat1 * matrix ;
            } else {
              // matrix = kron( Ileft , mat1 , Iright ) * matrix
              const int64_t nLeft  = 1 << qubit ;
              const int64_t nRight = 1 << ( nbQubits - qubit - 1 ) ;
              #pragma omp parallel for
              for ( int64_t j = 0; j < matrix.cols(); j++ ) {
                for ( int64_t k = 0; k < nLeft; k++ ) {
                  for ( int64_t i = 0; i < nRight; i++ ) {
                    const int64_t i1 = i + 2 * k * nRight ;
                    const int64_t i2 = i1 + nRight ;
                    const T x1 = matrix(i1,j) ;
                    const T x2 = matrix(i2,j) ;
                    matrix(i1,j) = mat1(0,0) * x1 + mat1(0,1) * x2 ;
                    matrix(i2,j) = mat1(1,0) * x1 + mat1(1,1) * x2 ;
                  }
                }
              }
            }
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

