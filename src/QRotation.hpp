//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_QRotation_hpp
#define qclab_QRotation_hpp

#include "QAngle.hpp"
#include "QAdjustable.hpp"
#include <cassert>

namespace qclab {

  /**
   * \class QRotation
   * \brief Class for representing an adjustable quantum rotation.
   *
   * This class stores the quantum angle \f$\theta/2\f$, i.e., the cosine and
   * sine of half the adjustable quantum rotation parameter \f$\theta\f$.
   */
  template <typename T>
  class QRotation : public QAdjustable
  {

    public:
      /// Quantum angle type of this quantum rotation.
      using angle_type = QAngle< T > ;

      /**
       * \brief Default constructor. Constructs an adjustable quantum rotation
       *        with angle \f$\theta = 0\f$.
       */
      QRotation()
      : angle_( 0 )
      , QAdjustable( false )
      { } // QRotation()

      /**
       * \brief Constructs an adjustable quantum rotation with the given
       *        quantum angle `angle` = \f$\theta/2\f$ and the flag `fixed`.
       *        The default value of `fixed` is false.
       */
      QRotation( const angle_type& angle , const bool fixed = false )
      : angle_( angle )
      , QAdjustable( fixed )
      { } // QRotation(angle,fixed)

      /**
       * \brief Constructs an adjustable quantum rotation with the given
       *        value `theta` = \f$\theta\f$ and the flag `fixed`.
       *        The default value of `fixed` is false.
       */
      QRotation( const T theta , const bool fixed = false )
      : angle_( theta/2 )
      , QAdjustable( fixed )
      { } // QRotation(theta,fixed)

      /**
       * \brief Constructs an adjustable quantum rotation with the given values
       *        `cos` = \f$\cos(\theta/2)\f$ and `sin` = \f$\sin(\theta/2)\f$
       *        and the flag `fixed`. The default value of `fixed` is false.
       */
      QRotation( const T cos , const T sin , const bool fixed = false )
      : angle_( cos , sin )
      , QAdjustable( fixed )
      { } // QRotation(cos,sin,fixed)

      /// Checks if this quantum rotation is fixed.
      inline bool fixed() const { return QAdjustable::fixed() ; }

      /// Returns the quantum angle \f$\theta/2\f$ of this quantum rotation.
      inline const angle_type& angle() const { return angle_ ; }

      /// Returns the numerical value \f$\theta\f$ of this quantum rotation.
      inline T theta() const { return 2 * angle_.theta() ; }

      /// Returns the cosine \f$\cos(\theta/2)\f$ of this quantum rotation.
      inline T cos() const { return angle_.cos() ; }

      /// Returns the sine \f$\sin(\theta/2)\f$ of this quantum rotation.
      inline T sin() const { return angle_.sin() ; }

      /**
       * \brief Updates this quantum rotation with the given quantum angle
       *        `angle` = \f$\theta/2\f$.
       */
      void update( const angle_type& angle ) {
        assert( !fixed() ) ;
        angle_.update( angle ) ;
      }

      /**
       * \brief Updates this quantum rotation with the given value
       *        `theta` = \f$\theta\f$.
       */
      void update( const T theta ) {
        assert( !fixed() ) ;
        angle_.update( theta/2 ) ;
      }

      /**
       * \brief Updates this quantum rotation with the given values
       *        `cos` = \f$\cos(\theta/2)\f$ and `sin` = \f$\sin(\theta/2)\f$.
       */
      void update( const T cos , const T sin ) {
        assert( !fixed() ) ;
        angle_.update( cos , sin ) ;
      }

      /// Checks if `other` is equal to this quantum rotation.
      inline bool operator==( const QRotation< T >& other ) const {
        return other.angle() == angle_ ;
      }

      /// Checks if `other` is different from this quantum rotation.
      inline bool operator!=( const QRotation< T >& other ) const {
        return other.angle() != angle_ ;
      }

    protected:
      angle_type  angle_ ;  ///< Quantum angle of this quantum rotation.

  } ; // class QRotation

} // namespace qclab

#endif

