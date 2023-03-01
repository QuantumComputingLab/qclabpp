//  (C) Copyright Roel Van Beeumen and Daan Camps 2021.

#pragma once

#include "qclab/qgates/CX.hpp"

namespace qclab {

  namespace qgates {

    template <typename T>
    using CNOT = CX< T > ;

  } // namespace qgates

} // namespace qclab

