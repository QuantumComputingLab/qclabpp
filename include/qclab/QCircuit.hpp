//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#pragma once

#include "qclab/QObject.hpp"
#include "qclab/QAdjustable.hpp"
#include "qclab/dense/kron.hpp"
#include "qclab/dense/transpose.hpp"
#include "qclab/io/QASMFile.hpp"
#include <cassert>
#include <numeric>
#include <vector>

namespace qclab {

  /**
   * \class QCircuit
   * \brief Class for representing a quantum circuit.
   *
   * This class stores the gates and sub-circuits of a quantum circuit.
   */
  template <typename T, typename G = QObject< T >>
  class QCircuit : public QObject< T > , public QAdjustable
  {

    public:
      /// Gate type of this quantum circuit.
      using gate_type = G ;
      /// Vector type of this quantum circuit.
      using vector_type = std::vector< std::unique_ptr< gate_type > > ;
      /// Size type of this quantum circuit.
      using size_type        = typename vector_type::size_type ;
      /// Reference type of this quantum circuit.
      using reference        = typename vector_type::reference ;
      /// Const reference type of this quantum circuit.
      using const_reference  = typename vector_type::const_reference ;
      /// Iterator type of this quantum circuit.
      using iterator         = typename vector_type::iterator ;
      /// Const iterator type of this quantum circuit.
      using const_iterator   = typename vector_type::const_iterator ;
      /// Reverse iterator type of this quantum circuit.
      using reverse_iterator = typename vector_type::reverse_iterator ;
      /// Const reverse iterator type of this quantum circuit.
      using const_reverse_iterator =
                               typename vector_type::const_reverse_iterator ;

      /// Constructs a quantum circuit of `nbQubits` qubits starting at qubit 0.
      QCircuit( const int nbQubits )
      : nbQubits_( nbQubits )
      , offset_( 0 )
      {
        assert( nbQubits > 0 ) ;
      } // QCircuit(nbQubits)

      /**
       * \brief Constructs a quantum circuit of `nbQubits` qubits starting at
       *        qubit `offset`.
       */
      QCircuit( const int nbQubits , const int offset )
      : nbQubits_( nbQubits )
      , offset_( offset )
      {
        assert( nbQubits > 0 ) ;
        assert( offset >= 0 ) ;
      } // QCircuit(nbQubits,offset)

      /**
       * \brief Constructs a quantum circuit of `nbQubits` qubits, consisting of
       *        `nbGates` gates and starting at qubit `offset`.
       */
      QCircuit( const int nbQubits , const int offset , const int nbGates )
      : nbQubits_( nbQubits )
      , offset_( offset )
      , gates_( nbGates )
      {
        assert( nbQubits > 0 ) ;
        assert( offset >= 0 ) ;
        assert( nbGates >= 0 ) ;
      } // QCircuit(nbQubits,offset,gates)

      /// Constructs a quantum circuit from the QASM file `filename`.
      QCircuit( const std::string filename , const int offset = 0 )
      : offset_( offset )
      {
        // QASM file
        qclab::io::QASMFile qasm( filename ) ;
        // number of qubits
        nbQubits_ = qasm.nbQubits() ;
        // loop over gates in file
        for ( auto& gate : qasm.gates< G >() ) {
          push_back( std::move( gate ) ) ;
        }
      } // QCircuit(filename)

      // nbQubits
      inline int nbQubits() const override { return nbQubits_ ; }

      // fixed
      inline bool fixed() const override { return QAdjustable::fixed() ; }

      // controlled
      inline bool controlled() const override { return false ; }

      // qubit
      inline int qubit() const override { return offset_ ; }

      // setQubit
      inline void setQubit( const int qubit ) override { assert( false ) ; }

      // qubits
      std::vector< int > qubits() const override {
        std::vector< int > v( nbQubits_ ) ;
        std::iota( v.begin() , v.end() , offset_ ) ;
        return v ;
      }

      // setQubits
      inline void setQubits( const int* qubits ) override { assert( false ) ; }

      // matrix
      qclab::dense::SquareMatrix< T > matrix() const override {
        auto mat = qclab::dense::eye< T >( 1 << nbQubits_ ) ;
        for ( auto it = begin(); it != end(); ++it ) {
          (*it)->apply( Side::Right , Op::NoTrans , nbQubits_ , mat ) ;
        }
        return mat ;
      }

      // apply
      void apply( Op op , const int nbQubits , std::vector< T >& vector ,
                  const int offset = 0 ) const override {
        if ( op == Op::NoTrans ) {
          // NoTrans
          for ( auto it = begin(); it != end(); ++it ) {
            (*it)->apply( op , nbQubits , vector , offset_ + offset ) ;
          }
        } else {
          // [Conj]Trans
          for ( auto it = rbegin(); it != rend(); ++it ) {
            (*it)->apply( op , nbQubits , vector , offset_ + offset ) ;
          }
        }
      }

    #ifdef QCLAB_OMP_OFFLOADING
      // apply_device
      void apply_device( Op op , const int nbQubits , T* vector ,
                         const int offset = 0 ) const override {
        if ( op == Op::NoTrans ) {
          // NoTrans
          for ( auto it = begin(); it != end(); ++it ) {
            (*it)->apply_device( op , nbQubits , vector , offset_ + offset ) ;
          }
        } else {
          // [Conj]Trans
          for ( auto it = rbegin(); it != rend(); ++it ) {
            (*it)->apply_device( op , nbQubits , vector , offset_ + offset ) ;
          }
        }
      }
    #endif

      /// Simulates this quantum circuit for the given vector `vector`.
      void simulate( std::vector< T >& vector ) const {
        apply( Op::NoTrans , nbQubits_ , vector ) ;
      }

    #ifdef QCLAB_OMP_OFFLOADING
      /// Simulates this quantum circuit for the given vector `vector`.
      void simulate_device( T* vector ) const {
        apply_device( Op::NoTrans , nbQubits_ , vector ) ;
      }
    #endif

      // apply
      void apply( Side side , Op op , const int nbQubits ,
                  qclab::dense::SquareMatrix< T >& matrix ,
                  const int offset = 0 ) const override {
        if ( op == Op::NoTrans ) {
          if ( side == Side::Left ) {
            // Left + NoTrans
            for ( auto it = rbegin(); it != rend(); ++it ) {
              (*it)->apply( side , op , nbQubits , matrix , offset_ + offset ) ;
            }
          } else {
            // Right + NoTrans
            for ( auto it = begin(); it != end(); ++it ) {
              (*it)->apply( side , op , nbQubits , matrix , offset_ + offset ) ;
            }
          }
        } else {
          if ( side == Side::Left ) {
            // Left + [Conj]Trans
            for ( auto it = begin(); it != end(); ++it ) {
              (*it)->apply( side , op , nbQubits , matrix , offset_ + offset ) ;
            }
          } else {
            // Right + [Conj]Trans
            for ( auto it = rbegin(); it != rend(); ++it ) {
              (*it)->apply( side , op , nbQubits , matrix , offset_ + offset ) ;
            }
          }
        }
      }

      // print
      void print() const override {
        printMatrix( matrix() ) ;
      }

      // toQASM
      int toQASM( std::ostream& stream , const int offset = 0 ) const override {
        for ( auto it = begin(); it != end(); ++it ) {
          int out = (*it)->toQASM( stream , offset_ + offset ) ;
          if ( out != 0 ) return out ;
        }
        return 0 ;
      }

      // operator==

      // operator!=

      // equals
      inline bool equals( const QObject< T >& other ) const override {
        return ( other.matrix() == matrix() ) ;
      }

      /// Returns the qubit offset of this quantum circuit.
      int offset() const { return offset_ ; }

      /// Sets the qubit offset of this quantum circuit.
      void setOffset( const int offset ) { offset_ = offset ; }

      //
      // Element access
      //

      /// Access specified gate of this quantum circuit with bounds checking.
      reference at( const size_type pos ) {
        assert( pos < gates_.size() ) ;
        return gates_.at( pos ) ;
      }

      /// Access specified gate of this quantum circuit with bounds checking.
      const_reference at( const size_type pos ) const {
        assert( pos < gates_.size() ) ;
        return gates_.at( pos ) ;
      }

      /// Access specified gate of this quantum circuit
      reference operator[]( const size_type pos ) {
        assert( pos < gates_.size() ) ;
        return gates_[ pos ] ;
      }

      /// Access specified gate of this quantum circuit
      const_reference operator[]( const size_type pos ) const {
        assert( pos < gates_.size() ) ;
        return gates_[ pos ] ;
      }

      /// Access the first gate of this quantum circuit.
      reference front() { return gates_.front() ; }

      /// Access the first gate of this quantum circuit.
      const_reference front() const { return gates_.front()() ; }

      /// Access the last gate of this quantum circuit.
      reference back() { return gates_.back() ; }

      /// Access the last gate of this quantum circuit.
      const_reference back() const { return gates_.back() ; }

      //
      // Iterators
      //

      /// Returns an iterator to the beginning of this quantum circuit.
      iterator begin() { return gates_.begin() ; }

      /// Returns a const iterator to the beginning of this quantum circuit.
      const_iterator begin() const { return gates_.begin() ; }

      /// Returns an iterator to the end of this quantum circuit.
      iterator end() { return gates_.end() ; }

      /// Returns a const iterator to the end of this quantum circuit.
      const_iterator end() const { return gates_.end() ; }

      /// Returns a reverse iterator to the beginning of this quantum circuit.
      reverse_iterator rbegin() { return gates_.rbegin() ; }

      /// Returns a const reverse iterator to the beginning of this quantum circuit.
      const_reverse_iterator rbegin() const { return gates_.rbegin() ; }

      /// Returns a reverse iterator to the end of this quantum circuit.
      reverse_iterator rend() { return gates_.rend() ; }

      /// Returns a const reverse iterator to the end of this quantum circuit.
      const_reverse_iterator rend() const { return gates_.rend() ; }

      //
      // Capacity
      //

      /// Checks if this quantum circuit is empty.
      bool empty() const { return gates_.empty() ; }

      /// Returns the number of gates in this quantum circuit.
      size_type nbGates() const { return gates_.size() ; }

      /// Returns the maximum possible number of gates in this quantum circuit.
      size_type maxNbGates() const { return gates_.max_size() ; }

      /// Reserves storage for this quantum circuit.
      void reserve( const size_type new_cap ) { gates_.reserve( new_cap ) ; }

      /**
       * \brief Returns the number of elements that can be held in currently
       *        allocated storage of this quantum circuit.
       */
      size_type capacity() const { return gates_.capacity() ; }

      /// Reduces memory usage by freeing unused memory of this quantum circuit.
      void shrink_to_fit() { gates_.shrink_to_fit() ; }

      //
      // Modifiers
      //

      /// Clears the gates of this quantum circuit.
      void clear() { gates_.clear() ; }

      /// Inserts gate `gate` at position `pos` in this quantum circuit.
      iterator insert( const_iterator pos , std::unique_ptr< gate_type > gate ){
        assert( canInsert( *gate ) ) ;
        return gates_.insert( pos , std::move( gate ) ) ;
      }

      /// Inserts gates at position `pos` in this quantum circuit.
      template <typename InputIt>
      iterator insert( const_iterator pos , InputIt first , InputIt last ) {
        if ( first == last ) {
          // return (non-const) iterator equal to pos
          iterator it = this->begin() ;
          std::advance( it , std::distance< const_iterator >( it , pos ) ) ;
          return it ;
        }
        iterator it = this->insert( pos , std::move( *first ) ) ;
        for ( auto itin = first + 1; itin != last; ++itin ) {
          it = this->insert( ++it , std::move( *itin ) ) ;
        }
        // return (non-const) iterator pointing to the first gate inserted
        std::advance( it , -std::distance( first + 1 , last ) ) ;
        return it ;
      }

      /// Erases the gate at position `pos` from this quantum circuit.
      iterator erase( const_iterator pos ) { return gates_.erase( pos ) ; }

      /**
       * \brief Erases the gates in the range [`first`,`last`) from this
       *        quantum circuit.
       */
      iterator erase( const_iterator first , const_iterator last ) {
        return gates_.erase( first , last ) ;
      }

      /// Adds a gate to the end of this quantum circuit.
      void push_back( std::unique_ptr< gate_type > gate ) {
        assert( canInsert( *gate ) ) ;
        gates_.push_back( std::move( gate ) ) ;
      }

      /// Removes the last gate of this quantum circuit.
      void pop_back() { gates_.pop_back() ; }

      /// Changes the number of gates of this quantum circuit stored.
      void resize( const size_type count ) { gates_.resize( count ) ; }

      /// Checks if the gate `gate` can be inserted into this quantum circuit.
      bool canInsert( const gate_type& gate ) const {
        const auto qubits = gate.qubits() ;
        if ( *std::max_element( qubits.begin() , qubits.end() ) >= nbQubits_ ) {
          return false ;
        }
        return true ;
      }

    protected:
      /// Number of qubits of this quantum circuit.
      int          nbQubits_ ;
      /// Qubit offset of this quantum circuit.
      int          offset_ ;
      /// Unique pointers to the gates of this quantum circuit.
      vector_type  gates_ ;

  } ; // class QCircuit


  /**
   * \brief Returns the Frobenius norm of the difference between the given
   *        circuits `circuit1` and `circuit2`.
   */
  template <typename C1, typename C2>
  double nrmF( const C1& circuit1 , const C2& circuit2 ) {
    assert( circuit1.nbQubits() == circuit2.nbQubits() ) ;
    const auto m1 = circuit1.matrix() ;
    const auto m2 = circuit2.matrix() ;
    double scale = 0.0 ;
    #pragma omp parallel for reduction(max:scale)
    for ( int c = 0; c < m1.cols(); c++ ) {
      for ( int r = 0; r < m1.rows(); r++ ) {
        const double diff = std::abs(  m1(r,c) - m2(r,c) ) ;
        scale = std::max( scale , diff ) ;
      }
    }
    if ( scale == 0.0 ) { scale = 1.0 ; }
    double trace = 0.0 ;
    #pragma omp parallel for reduction(+:trace)
    for ( int c = 0; c < m1.cols(); c++ ) {
      for ( int r = 0; r < m1.rows(); r++ ) {
        const double real = std::real( m1(r,c) - m2(r,c) ) / scale ;
        const double imag = std::imag( m1(r,c) - m2(r,c) ) / scale ;
        trace += real * real + imag * imag ;
      }
    }
    return scale * std::sqrt( trace ) ;
  }

} // namespace qclab

