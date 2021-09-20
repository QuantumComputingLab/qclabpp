cmake_minimum_required( VERSION 3.16 FATAL_ERROR )

# project
project( QCLABPP VERSION 0.1.2 LANGUAGES CXX )

# fetch content
include( FetchContent )

# GTest
FetchContent_Declare( gtest
  GIT_REPOSITORY https://github.com/google/googletest.git
)
FetchContent_MakeAvailable( gtest )

# qclab/include
add_library( qclabpp INTERFACE )
target_include_directories( qclabpp INTERFACE ${PROJECT_SOURCE_DIR}/include )
target_compile_features( qclabpp INTERFACE cxx_std_17 )

# qclab/test
add_subdirectory( test )

# qclab/examples
add_subdirectory( examples )

# openmp
find_package( OpenMP )
if( TARGET OpenMP::OpenMP_CXX )
  target_link_libraries( qclabpp INTERFACE OpenMP::OpenMP_CXX )
endif()

# documentation
find_package( Doxygen )
if( DOXYGEN_FOUND )
  configure_file( ${CMAKE_SOURCE_DIR}/doc/doxygen/doxygen.dox.in
    ${CMAKE_BINARY_DIR}/doxygen.dox @ONLY)
  add_custom_target( doc ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/doxygen.dox
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating API documentation with doxygen" VERBATIM )
endif()

