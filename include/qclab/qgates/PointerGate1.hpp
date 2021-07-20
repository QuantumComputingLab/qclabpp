//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_PointerGate1_hpp
#define qclab_qgates_PointerGate1_hpp

#include "qclab/qgates/QGate1.hpp"

namespace qclab {

  namespace qgates {

    /**
     * \class PointerGate1
     * \brief 1-qubit gate pointing to another 1-qubit gate.
     */
    template <typename T>
    class PointerGate1 : public QGate1< T >
    {

      public:
        /**
         * \brief Constructs a 1-qubit pointer gate from the given gate pointer
         *        `gate` and qubit offset `offset`.
         *        The default value of `offset` is 0.
         */
        PointerGate1( QGate1< T >* gate , const int offset = 0 )
        : offset_( offset )
        , gate_( gate )
        { } // PointerGate1(gate,offset)

        // nbQubits

        // fixed
        inline bool fixed() const override { return gate_->fixed() ; }

        // controlled

        // qubit
        inline int qubit() const override { return gate_->qubit() + offset_ ; }

        // setQubit
        inline void setQubit( const int qubit ) override { assert( false ) ; }

        // qubits
        std::vector< int > qubits() const override {
          return std::vector< int >( { qubit() } ) ;
        }

        // setQubits
        inline void setQubits( const int* qubits ) override { assert( false ) ;}

        // matrix
        qclab::dense::SquareMatrix< T > matrix() const override {
          return gate_->matrix() ;
        }

        // apply

        // print

        // toQASM
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          return gate_->toQASM( stream , offset_ + offset ) ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using G = PointerGate1< T > ;
          if ( const G* p = dynamic_cast< const G* >( &other ) ) {
            return *p->ptr() == *this->ptr() ;
          }
          return other == *this->ptr() ;
        }

        /// Returns the qubit offset of this 1-qubit pointer gate.
        int offset() const { return offset_ ; }

        /// Sets the qubit offset of this 1-qubit pointer gate.
        void setOffset( const int offset ) { offset_ = offset ; }

        /// Returns a pointer to the gate of this 1-qubit pointer gate.
        QGate1< T >* ptr() { return gate_ ; }

        /// Returns a const pointer to the gate of this 1-qubit pointer gate.
        const QGate1< T >* ptr() const { return gate_ ; }

      protected:
        int           offset_ ;  ///< Qubit offset of this 1-qubit pointer gate.
        QGate1< T >*  gate_ ;    ///< Gate pointer of this 1-qubit pointer gate.

    } ; // class PointerGate1

  } // namespace qgates

} // namespace qclab

#endif

