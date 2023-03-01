#include <gtest/gtest.h>

int main( int argc , char **argv ) {

#ifdef QCLAB_OMP_OFFLOADING
  std::cout << "\n!!! QCLAB OMP OFFLOADING !!!\n" << std::endl ;
#endif

  ::testing::InitGoogleTest( &argc , argv ) ;
  return RUN_ALL_TESTS() ;

}

