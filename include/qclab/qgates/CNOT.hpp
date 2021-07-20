//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#ifndef qclab_qgates_CNOT_hpp
#define qclab_qgates_CNOT_hpp

#include "qclab/qgates/CX.hpp"

namespace qclab {

  namespace qgates {

    template <typename T>
    using CNOT = CX< T > ;

  } // namespace qgates

} // namespace qclab

#endif

