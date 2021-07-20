//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_QRotation_hpp
#define qclab_QRotation_hpp

#include "qclab/QAngle.hpp"

namespace qclab {

  /**
   * \class QRotation
   * \brief Class for representing a quantum rotation.
   *
   * This class stores the quantum angle \f$\theta/2\f$, i.e., the cosine and
   * sine of half the quantum rotation parameter \f$\theta\f$.
   */
  template <typename T>
  class QRotation
  {

    public:
      /// Quantum angle type of this quantum rotation.
      using angle_type = QAngle< T > ;

      /**
       * \brief Default constructor. Constructs a aquantum rotation with angle
       *        \f$\theta = 0\f$.
       */
      QRotation()
      : angle_( 0 )
      { } // QRotation()

      /**
       * \brief Constructs a quantum rotation with the given quantum angle
       *        `angle` = \f$\theta/2\f$.
       */
      QRotation( const angle_type& angle )
      : angle_( angle )
      { } // QRotation(angle)

      /**
       * \brief Constructs a quantum rotation with the given value
       *        `theta` = \f$\theta\f$.
       */
      QRotation( const T theta )
      : angle_( theta/2 )
      { } // QRotation(theta)

      /**
       * \brief Constructs a quantum rotation with the given values
       *        `cos` = \f$\cos(\theta/2)\f$ and `sin` = \f$\sin(\theta/2)\f$.
       */
      QRotation( const T cos , const T sin )
      : angle_( cos , sin )
      { } // QRotation(cos,sin)

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
        angle_.update( angle ) ;
      }

      /**
       * \brief Updates this quantum rotation with the given value
       *        `theta` = \f$\theta\f$.
       */
      void update( const T theta ) {
        angle_.update( theta/2 ) ;
      }

      /**
       * \brief Updates this quantum rotation with the given values
       *        `cos` = \f$\cos(\theta/2)\f$ and `sin` = \f$\sin(\theta/2)\f$.
       */
      void update( const T cos , const T sin ) {
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

      /// Multiplies `rhs` to this quantum rotation.
      inline QRotation< T >& operator*=( const QRotation< T >& rhs ) {
        update( angle_ + rhs.angle() ) ;
        return *this ;
      }

      /// Multiplies the inverse of `rhs` to this quantum rotation.
      inline QRotation< T >& operator/=( const QRotation< T >& rhs ) {
        update( angle_ - rhs.angle() ) ;
        return *this ;
      }

      /// Multiplies 2 quantum rotations.
      friend QRotation< T > operator*( QRotation< T > lhs ,
                                       const QRotation< T >& rhs ) {
        lhs *= rhs ;
        return lhs ;
      }

      /// Multiplies `lhs` and the inverse of `rhs`.
      friend QRotation< T > operator/( QRotation< T > lhs ,
                                       const QRotation< T >& rhs ) {
        lhs /= rhs ;
        return lhs ;
      }

      /// Returns the inverse this quantum rotation.
      inline QRotation< T > inv() const {
        QRotation< T > rotation( -angle_ ) ;
        return rotation ;
      }

    protected:
      angle_type  angle_ ;  ///< Quantum angle of this quantum rotation.

  } ; // class QRotation

} // namespace qclab

#endif

