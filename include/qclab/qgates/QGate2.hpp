//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#pragma once

#include "qclab/QObject.hpp"
#include "qclab/dense/transpose.hpp"

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
          printMatrix4x4( this->matrix() ) ;
        }

        // toQASM

        // operator==

        // operator!=

        // equals

    } ; // class QGate2

  } // namespace qgates

} // namespace qclab

