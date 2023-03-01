//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#pragma once

#include "qclab/QObject.hpp"
#include "qclab/dense/transpose.hpp"

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

