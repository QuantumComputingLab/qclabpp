//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_QObject_hpp
#define qclab_QObject_hpp

#include "util.hpp"
#include "qasm.hpp"
#include "dense/SquareMatrix.hpp"
#include <vector>

/**
 * QCLAB++ namespace.
 */
namespace qclab {

  /**
   * \class QObject
   * \brief General base class for all quantum objects.
   */
  template <typename T>
  class QObject
  {

    public:
      /// Value type of this quantum object.
      using value_type = T ;

      /// Returns the number of qubits of this quantum object.
      virtual int nbQubits() const = 0 ;

      /// Checks if this quantum object is fixed.
      virtual bool fixed() const = 0 ;

      /// Checks if this quantum object is controlled.
      virtual bool controlled() const = 0 ;

      /// Returns the first qubit of this quantum object.
      virtual int qubit() const = 0 ;

      /// Sets the qubit of 1-qubit quantum objects, otherwise void.
      virtual void setQubit( const int qubit ) = 0 ;

      /// Returns the qubits of this quantum object in ascending order.
      virtual std::vector< int > qubits() const = 0 ;

      /// Sets the qubits of this quantum object.
      virtual void setQubits( const int* qubits ) = 0 ;

      /// Returns the unitary matrix corresponding to this quantum object.
      virtual qclab::dense::SquareMatrix< T > matrix() const = 0 ;

      /// Applies this quantum object to the given matrix.
      virtual void apply( Side side , Op op , const int size ,
                          qclab::dense::SquareMatrix< T >& matrix ,
                          const int offset = 0 ) const = 0 ;

      /// Prints this quantum object.
      virtual void print() const = 0 ;

      /// Writes the QASM code of this quantum object to the given `stream`.
      virtual int toQASM( std::ostream& stream ,
                          const int offset = 0 ) const = 0 ;

      /// Checks if `other` is equal to this quantum object.
      inline bool operator==( const QObject< T >& other ) const {
        return equals( other ) ;
      }

      /// Checks if `other` is different from this quantum object.
      inline bool operator!=( const QObject< T >& other ) const {
        return !( *this == other ) ;
      }

    protected:
      /// Checks if `other` equals this quantum object.
      virtual bool equals( const QObject< T >& other ) const = 0 ;

  } ; // class QObject

} // namespace qclab

#endif

