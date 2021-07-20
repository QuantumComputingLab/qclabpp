//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_Phase45_hpp
#define qclab_qgates_Phase45_hpp

#include "qclab/qgates/QGate1.hpp"

namespace qclab {

  namespace qgates {

    /**
     * \class Phase45
     * \brief 1-qubit T gate (pi/4 rotation about the Z-axis).
     */
    template <typename T>
    class Phase45 : public QGate1< T >
    {

      public:
        /// Default constructor. Constructs a T gate on qubit 0.
        Phase45()
        : QGate1< T >()
        { } // Phase45()

        /// Constructs a T gate on the given qubit `qubit`.
        Phase45( const int qubit )
        : QGate1< T >( qubit )
        { } // Phase45(qubit)

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
          using R = qclab::real_t< T > ;
          const R sqrt2 = R(1) / std::sqrt( R(2) ) ;
          return qclab::dense::SquareMatrix< T >( 1 ,       0        ,
                                                  0 , T(sqrt2,sqrt2) ) ;
        }

        // apply

        // print

        // toQASM
        int toQASM( std::ostream& stream ,
                    const int offset = 0 ) const override {
          stream << qasmT( this->qubit_ + offset ) ;
          return 0 ;
        }

        // operator==

        // operator!=

        // equals
        inline bool equals( const QObject< T >& other ) const override {
          using G = Phase45< T > ;
          if ( const G* p = dynamic_cast< const G* >( &other ) ) return true ;
          return false ;
        }

    } ; // class Phase45

  } // namespace qgates

} // namespace qclab

#endif

