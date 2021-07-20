//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_Identity_hpp
#define qclab_qgates_Identity_hpp

#include "qclab/qgates/QGate1.hpp"

namespace qclab {

  namespace qgates {

    /**
     * \class Identity
     * \brief 1-qubit Identity gate.
     */
    template <typename T>
    class Identity : public QGate1< T >
    {

      public:
        /// Default constructor. Constructs an identity gate on qubit 0.
        Identity()
        : QGate1< T >()
        { } // Identity()

        /// Constructs an identity gate on the given qubit `qubit`.
        Identity( const int qubit )
        : QGate1< T >( qubit )
        { } // Identity(qubit)

        // nbQubits

        // fixed
        inline bool fixed() const override { return true ; } ;

        // controlled

        // qubit

        // setQubit

        // qubits

        // setQubits

        // matrix
        qclab::dense::SquareMatrix< T > matrix() const override {
          return qclab::dense::SquareMatrix< T >( 1 , 0 ,
                                                  0 , 1 ) ;
        }

        // apply

        // print

        // toQASM
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          // do nothing
          return 0 ;
        }

        // operator==

        // operator!=

      protected:
        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using I = Identity< T > ;
          if ( const I* p = dynamic_cast< const I* >( &other ) ) return true ;
          return false ;
        }

    } ; // class Identity

  } // namespace qgates

} // namespace qclab

#endif

