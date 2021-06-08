cmake_minimum_required( VERSION 3.16 FATAL_ERROR )

# project
project( QCLABPP VERSION 0.1.1 LANGUAGES CXX )
set( CMAKE_CXX_STANDARD 14 )

# openmp
find_package( OpenMP )

# fetch content
include( FetchContent )

# GTest
FetchContent_Declare( gtest
  GIT_REPOSITORY https://github.com/google/googletest.git
)
FetchContent_MakeAvailable( gtest )

# qclab/src
add_library( qclab INTERFACE IMPORTED )
target_include_directories( qclab INTERFACE ${CMAKE_SOURCE_DIR}/src )

# qclab/test
add_subdirectory( test )

# qclab/examples
add_subdirectory( examples )

# documentation
find_package( Doxygen )
if( DOXYGEN_FOUND )
  configure_file( ${CMAKE_SOURCE_DIR}/doc/doxygen/doxygen.dox.in
    ${CMAKE_BINARY_DIR}/doxygen.dox @ONLY)
  add_custom_target( doc ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/doxygen.dox
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating API documentation with doxygen" VERBATIM )
endif()

