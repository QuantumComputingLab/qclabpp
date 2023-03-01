//  (C) Copyright Roel Van Beeumen 2021.

#pragma once

#include "qclab/dense/memory.hpp"

namespace qclab {

  /**
   * Namespace qclab::dense.
   */
  namespace dense {

    /**
     * \class SquareMatrix
     * \brief Class for representing a square matrix.
     *
     * This class stores the size and data of a square matrix.
     */
    template <typename T>
    class SquareMatrix
    {

      public:
        /// Value type of this square matrix.
        using value_type = T ;
        /// Size type of this square matrix.
        using size_type  = int64_t ;
        /// Data type of this square matrix.
        using data_type  = std::unique_ptr< T[] > ;
        /// Pointer type of this square matrix.
        using ptr_type   = T* ;

        /// Default constructor. Constructs a square matrix of size 1.
        SquareMatrix()
        : size_( 1 )
        , data_( alloc_unique_array< T >( 1 ) )
        { } // SquareMatrix()

        /// Constructs a square matrix of size `size`.
        SquareMatrix( const int64_t size )
        : size_( size )
        , data_( alloc_unique_array< T >( size * size ) )
        {
          assert( size > 0 ) ;
        } // SquareMatrix(size)

        /// Constructs a square matrix of size `size` initialize with `value`.
        SquareMatrix( const int64_t size , const T value )
        : size_( size )
        , data_( init_unique_array< T >( size * size , value ) )
        {
          assert( size > 0 ) ;
        } // SquareMatrix(size,value)

        /// Constructs a 2x2 square matrix with the given elements.
        SquareMatrix( const T m00 , const T m01 ,
                      const T m10 , const T m11 )
        : size_( 2 )
        , data_( alloc_unique_array< T >( 4 ) )
        {
          data_[0] = m00 ; data_[1] = m10 ;
          data_[2] = m01 ; data_[3] = m11 ;
        } // SquareMatrix(m00,m10,m01,m11)

        /// Constructs a 4x4 square matrix with the given elements.
        SquareMatrix( const T m00 , const T m01 , const T m02 , const T m03 ,
                      const T m10 , const T m11 , const T m12 , const T m13 ,
                      const T m20 , const T m21 , const T m22 , const T m23 ,
                      const T m30 , const T m31 , const T m32 , const T m33 )
        : size_( 4 )
        , data_( alloc_unique_array< T >( 16 ) )
        {
          data_[ 0] = m00 ; data_[ 1] = m10 ; data_[ 2] = m20 ; data_[ 3] = m30;
          data_[ 4] = m01 ; data_[ 5] = m11 ; data_[ 6] = m21 ; data_[ 7] = m31;
          data_[ 8] = m02 ; data_[ 9] = m12 ; data_[10] = m22 ; data_[11] = m32;
          data_[12] = m03 ; data_[13] = m13 ; data_[14] = m23 ; data_[15] = m33;
        } // SquareMatrix(m00,m10,m20,m30,...,m33)

        /// Copy constructor.
        SquareMatrix( const SquareMatrix< T >& matrix )
        : size_( matrix.size() )
        , data_( alloc_unique_array< T >( size_ * size_ ) )
        {
          const T* data = matrix.ptr() ;
          #pragma omp parallel for
          for ( int64_t i = 0; i < size_*size_; i++ ) {
            data_[i] = data[i] ;
          }
        } // SquareMatrix(matrix)

        /// Copy assignment operator.
        SquareMatrix< T >& operator=( const SquareMatrix< T >& matrix ) {
          if ( matrix.size() != size_ ) {
            size_ = matrix.size() ;
            data_ = std::move( alloc_unique_array< T >( size_ * size_ ) ) ;
          }
          const T* data = matrix.ptr() ;
          #pragma omp parallel for
          for ( int64_t i = 0; i < size_*size_; i++ ) {
            data_[i] = data[i] ;
          }
          return *this ;
        }

        /// Returns the size of this square matrix.
        inline int64_t size() const { return size_ ; } ;

        /// Returns the number of rows of this square matrix.
        inline int64_t rows() const { return size_ ; } ;

        /// Returns the number of columns of this square matrix.
        inline int64_t cols() const { return size_ ; } ;

        /// Returns the leading dimension of this square matrix.
        inline int64_t ld() const { return size_ ; } ;

        /// Returns a pointer to this square matrix.
        inline T* ptr() { return data_.get() ; } ;

        /// Returns a const pointer to this square matrix.
        inline const T* ptr() const { return data_.get() ; } ;

        /// Returns the value of this square matrix at row `i` and column `j`.
        inline T& operator()( const int64_t i , const int64_t j ) {
          return data_[ i + j*size_ ] ;
        }

        /// Returns the value of this square matrix at row `i` and column `j`.
        inline const T& operator()( const int64_t i , const int64_t j ) const {
          return data_[ i + j*size_ ] ;
        }

        /// Checks if `other` is equal to this square matrix.
        inline bool operator==( const SquareMatrix< T >& other ) const {
          if ( other.size() != size_ ) return false ;
          for ( size_type j = 0; j < size_; j++ ) {
            for ( size_type i = 0; i < size_; i++ ) {
              if ( other(i,j) != (*this)(i,j) ) return false ;
            }
          }
          return true ;
        }

        /// Checks if `other` is different from this square matrix.
        inline bool operator!=( const SquareMatrix< T >& other ) const {
          return !( *this == other ) ;
        }

        /// Adds `rhs` to this square matrix.
        inline SquareMatrix< T >& operator+=( const SquareMatrix< T >& rhs ) {
          assert( rhs.size() == size_ ) ;
          #pragma omp parallel for
          for ( int64_t j = 0; j < size_; j++ ) {
            for ( int64_t i = 0; i < size_; i++ ) {
              (*this)(i,j) += rhs(i,j) ;
            }
          }
          return *this ;
        }

        /// Substracts `rhs` from this square matrix.
        inline SquareMatrix< T >& operator-=( const SquareMatrix< T >& rhs ) {
          assert( rhs.size() == size_ ) ;
          #pragma omp parallel for
          for ( int64_t j = 0; j < size_; j++ ) {
            for ( int64_t i = 0; i < size_; i++ ) {
              (*this)(i,j) -= rhs(i,j) ;
            }
          }
          return *this ;
        }

        /// Multiplies `rhs` to the right of this square matrix.
        inline SquareMatrix< T >& operator*=( const SquareMatrix< T >& rhs ) {
          assert( size_ == rhs.size() ) ;
          auto new_data = init_unique_array< T >( size_ * size_ ) ;
          #pragma omp parallel for
          for ( int64_t j = 0; j < size_; j++ ) {
            for ( int64_t i = 0; i < size_; i++ ) {
              for ( int64_t k = 0; k < size_; k++ ) {
                new_data.get()[ i + j*size_ ] += (*this)(i,k) * rhs(k,j) ;
              }
            }
          }
          data_ = std::move( new_data ) ;
          return *this ;
        }

        /// Adds 2 square matrices.
        friend SquareMatrix< T > operator+( SquareMatrix< T > lhs ,
                                            const SquareMatrix< T >& rhs ) {
          lhs += rhs ;
          return lhs ;
        }

        /// Substracts 2 square matrices.
        friend SquareMatrix< T > operator-( SquareMatrix< T > lhs ,
                                            const SquareMatrix< T >& rhs ) {
          lhs -= rhs ;
          return lhs ;
        }

        /// Multiplies 2 square matrices
        friend SquareMatrix< T > operator*( SquareMatrix< T > lhs ,
                                            const SquareMatrix< T >& rhs ) {
          lhs *= rhs ;
          return lhs ;
        }

      private:
        size_type  size_ ;
        data_type  data_ ;

    } ; // class SquareMatrix


    /// Returns a zero square matrix of size `size`.
    template <typename T>
    SquareMatrix< T > zeros( const int64_t size ) {
      return SquareMatrix< T >( size , T(0) ) ;
    }

    /// Returns an identity matrix of size `size`.
    template <typename T>
    SquareMatrix< T > eye( const int64_t size ) {
      SquareMatrix< T > mat( size , T(0) ) ;
      #pragma omp parallel for
      for ( int64_t i = 0; i < size; i++ ) {
        mat(i,i) = T(1) ;
      }
      return mat ;
    }

  } // namespace dense

} // namespace qclab

