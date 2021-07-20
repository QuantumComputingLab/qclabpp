//  (C) Copyright Roel Van Beeumen 2021.

#ifndef qclab_dense_memory_hpp
#define qclab_dense_memory_hpp

#include <cassert>
#include <memory>
#include <algorithm>

namespace qclab {

  namespace dense {

    /**
     * \brief Allocates an array of size `size`.
     */
    template <typename T>
    inline auto alloc_unique_array( const int64_t size ) {
      assert( size > 0 ) ;
      std::unique_ptr< T[] > mem( new T[size] ) ;
      assert( mem.get() != nullptr ) ;
      return mem ;
    } // alloc_unique_array

    /**
     * \brief Allocates an array of size `size` and intializes using the
     *        default constructor.
     */
    template <typename T>
    inline auto init_unique_array( const int64_t size ) {
      assert( size > 0 ) ;
      std::unique_ptr< T[] > mem( new T[size]() ) ;
      assert( mem.get() != nullptr ) ;
      return mem ;
    } // init_unique_array

    /**
     * \brief Allocates an array of size `size` and intializes with `value`.
     */
    template <typename T>
    inline auto init_unique_array( const int64_t size , const T value ) {
      auto mem = alloc_unique_array< T >( size ) ;
      std::fill( mem.get() , mem.get() + size , value ) ;
      return mem ;
    } // init_unique_array

  } // namespace qclab

} // namespace dense

#endif

