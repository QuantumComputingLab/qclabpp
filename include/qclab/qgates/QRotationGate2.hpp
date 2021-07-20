//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_QRotationGate2_hpp
#define qclab_qgates_QRotationGate2_hpp

#include "qclab/qgates/QGate2.hpp"
#include "qclab/QAdjustable.hpp"
#include "qclab/QRotation.hpp"
#include <array>

namespace qclab {

  namespace qgates {

    /**
     * \class QRotationGate2
     * \brief Base class for 2-qubit rotation gates.
     */
    template <typename T>
    class QRotationGate2 : public QGate2< T > , public QAdjustable
    {

      public:
        /// Real value type of this 2-qubit rotation gate.
        using real_type = qclab::real_t< T > ;
        /// Quantum angle type of this 2-qubit rotation gate.
        using angle_type = qclab::QAngle< real_type > ;

        /**
         * \brief Default constructor. Constructs a 2-qubit rotation gate on
         *        qubits 0 and 1 with parameter \f$\theta = 0\f$.
         */
        QRotationGate2()
        : qubits_( { 0 , 1 } )
        , rotation_()
        , QAdjustable()
        { } // QRotationGate2()

        /**
         * \brief Constructs a 2-qubit rotation gate on the given qubits
         *        `qubits` with quantum angle `angle` = \f$\theta/2\f$ and
         *        flag `fixed`. The default value of `fixed` is false.
         */
        QRotationGate2( const int* qubits , const angle_type& angle ,
                        const bool fixed = false )
        : rotation_( angle )
        , QAdjustable( fixed )
        {
          setQubits( qubits ) ;
        } // QRotationGate2(qubits,angle,fixed)

        /**
         * \brief Constructs a 2-qubit rotation gate on the given qubits
         *        `qubits` with value `theta` = \f$\theta\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        QRotationGate2( const int* qubits , const real_type theta ,
                        const bool fixed = false )
        : rotation_( theta )
        , QAdjustable( fixed )
        {
          setQubits( qubits ) ;
        } // QRotationGate2(qubits,theta,fixed)

        /**
         * \brief Constructs a 2-qubit rotation gate on the given qubits
         *        `qubits` with values `cos` = \f$\cos(\theta/2)\f$ and
         *        `sin` = \f$\sin(\theta/2)\f$, and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        QRotationGate2( const int* qubits , const real_type cos ,
                        const real_type sin , const bool fixed = false )
        : rotation_( cos , sin )
        , QAdjustable( fixed )
        {
          setQubits( qubits ) ;
        } // QRotationGate2(qubits,cos,sin,fixed)

        /**
         * \brief Constructs a 2-qubit rotation gate on the given
         *        qubits `qubit0` and `qubit1` with quantum angle
         *        `angle` = \f$\theta/2\f$ and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        QRotationGate2( const int qubit0 , const int qubit1 ,
                        const angle_type& angle , const bool fixed = false )
        : rotation_( angle )
        , QAdjustable( fixed )
        {
          const int qubits[] = { qubit0 , qubit1 } ;
          setQubits( &qubits[0] ) ;
        } // QRotationGate2(qubit0,qubit1,angle,fixed)

        /**
         * \brief Constructs a 2-qubit rotation gate on the given qubits
         *        `qubit0` and `qubit1` with value `theta` = \f$\theta\f$
         *        and flag `fixed`. The default value of `fixed` is false.
         */
        QRotationGate2( const int qubit0 , const int qubit1 ,
                        const real_type theta , const bool fixed = false )
        : rotation_( theta )
        , QAdjustable( fixed )
        {
          const int qubits[] = { qubit0 , qubit1 } ;
          setQubits( &qubits[0] ) ;
        } // QRotationGate2(qubit0,qubit1,theta,fixed)

        /**
         * \brief Constructs a 2-qubit rotation gate on the given qubits
         *        `qubit0` and `qubit1` with values `cos` = \f$\cos(\theta/2)\f$
         *        and `sin` = \f$\sin(\theta/2)\f$, and flag `fixed`.
         *        The default value of `fixed` is false.
         */
        QRotationGate2( const int qubit0 , const int qubit1 ,
                        const real_type cos , const real_type sin ,
                        const bool fixed = false )
        : rotation_( cos , sin )
        , QAdjustable( fixed )
        {
          const int qubits[] = { qubit0 , qubit1 } ;
          setQubits( &qubits[0] ) ;
        } // QRotationGate2(qubit0,qubit1,cos,sin,fixed)

        // nbQubits

        // fixed
        inline bool fixed() const override { return QAdjustable::fixed() ; }

        // controlled
        inline bool controlled() const override { return false ; }

        // qubit
        inline int qubit() const override { return qubits_[0] ; }

        // setQubit

        // qubits
        std::vector< int > qubits() const override {
          return std::vector< int >( { qubits_[0] , qubits_[1] } ) ;
        }

        // setQubits
        inline void setQubits( const int* qubits ) override {
          assert( qubits[0] >= 0 ) ; assert( qubits[1] >= 0 ) ;
          assert( qubits[0] != qubits[1] ) ;
          qubits_[0] = std::min( qubits[0] , qubits[1] ) ;
          qubits_[1] = std::max( qubits[0] , qubits[1] ) ;
        }

        // matrix

        // apply

        // print

        // toQASM

        // operator==

        // operator!=

        // equals

        /// Returns the quantum angle \f$\theta/2\f$ of this rotation gate.
        inline const angle_type& angle() const { return rotation_.angle() ; }

        /// Returns the numerical value \f$\theta\f$ of this rotation gate.
        inline real_type theta() const { return rotation_.theta() ; }

        /// Returns the cosine \f$\cos(\theta/2)\f$ of this rotation gate.
        inline real_type cos() const { return rotation_.cos() ; }

        /// Returns the sine \f$\sin(\theta/2)\f$ of this rotation gate.
        inline real_type sin() const { return rotation_.sin() ; }

        /**
         * \brief Updates this rotation gate with the given quantum angle
         *        `angle` = \f$\theta/2\f$.
         */
        void update( const angle_type& angle ) {
          assert( !fixed() ) ;
          rotation_.update( angle ) ;
        }

        /**
         * \brief Updates this rotation gate with the given value
         *        `theta` = \f$\theta\f$.
         */
        void update( const real_type theta ) {
          assert( !fixed() ) ;
          rotation_.update( theta ) ;
        }

        /**
         * \brief Updates this rotation gate with the given values
         *        `cos` = \f$\cos(\theta/2)\f$ and `sin` = \f$\sin(\theta/2)\f$.
         */
        void update( const real_type cos , const real_type sin ) {
          assert( !fixed() ) ;
          rotation_.update( cos , sin ) ;
        }

      protected:
        /// Qubits of this 2-qubit rotation gate.
        std::array< int , 2 >   qubits_ ;
        /// Quantum rotation of this 2-qubit rotation gate.
        QRotation< real_type >  rotation_ ;

    } ; // class QRotationGate2

  } // namespace qgates

} // namespace qclab

#endif

