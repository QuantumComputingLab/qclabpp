//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_QAngle_hpp
#define qclab_QAngle_hpp

#include <cassert>
#include <cmath>

namespace qclab {

  /**
   * \class QAngle
   * \brief Class for representing a quantum angle.
   *
   * This class stores the cosine and sine of the angle of a quantum angle.
   */
  template <typename T>
  class QAngle
  {

    public:
      /// Value type of this quantum angle.
      using value_type = T ;

      /// Default constructor. Constructs a quantum angle with \f$\theta = 0\f$.
      QAngle()
      : cos_( 1 )
      , sin_( 0 )
      { } // QAngle()

      /// Constructs a quantum angle with the given value `theta`.
      QAngle( const T theta )
      : cos_( std::cos( theta ) )
      , sin_( std::sin( theta ) )
      { } // QAngle(theta)

      /// Constructs a quantum angle with the given `cos` and `sin` values.
      QAngle( const T cos , const T sin )
      : cos_( cos )
      , sin_( sin )
      {
        const T eps = std::numeric_limits< T >::epsilon() ;
        assert( cos * cos + sin * sin - T(1) < 10*eps ) ;
      } // QAngle(cos,sin)

      /// Returns the numerical value \f$\theta\f$ of this quantum angle.
      inline const T theta() const { return std::atan2( sin_ , cos_ ) ; }

      /// Returns the cosine \f$\cos(\theta)\f$ of this quantum angle.
      inline const T cos() const { return cos_ ; }

      /// Returns the sine \f$\sin(\theta)\f$ of this quantum angle.
      inline const T sin() const { return sin_ ; }

      /// Updates this quantum angle with the given `angle`.
      void update( const QAngle< T >& angle ) {
        cos_ = angle.cos() ;
        sin_ = angle.sin() ;
      }

      /// Updates this quantum angle with the given value `theta`.
      void update( const T theta ) {
        cos_ = std::cos( theta ) ;
        sin_ = std::sin( theta ) ;
      }

      /// Updates this quantum angle with the given `cos` and `sin` values.
      void update( const T cos , const T sin ) {
        const T eps = std::numeric_limits< T >::epsilon() ;
        assert( cos * cos + sin * sin - T(1) < 10*eps ) ;
        cos_ = cos ;
        sin_ = sin ;
      }

      /// Checks if `other` is equal to this quantum angle.
      inline bool operator==( const QAngle< T >& other ) const {
        return ( ( other.cos() == cos_ ) && ( other.sin() == sin_ ) ) ;
      }

      /// Checks if `other` is different from this quantum angle.
      inline bool operator!=( const QAngle< T >& other ) const {
        return ( ( other.cos() != cos_ ) || ( other.sin() != sin_ ) ) ;
      }

    protected:
      T  cos_ ;  ///< Cosine value \f$\cos(\theta)\f$ of this quantum angle.
      T  sin_ ;  ///<   Sine value \f$\sin(\theta)\f$ of this quantum angle.

  } ; // class QAngle

} // namespace qclab

#endif

