add_library( qclabpp qgates/QGate1.cpp
                     qgates/Hadamard.cpp
                     qgates/Identity.cpp
                     qgates/PauliX.cpp
                     qgates/PauliY.cpp
                     qgates/PauliZ.cpp
                     qgates/RotationX.cpp
                     qgates/RotationY.cpp
                     qgates/RotationZ.cpp
                     qgates/Phase.cpp
                     qgates/Phase45.cpp
                     qgates/Phase90.cpp
                     qgates/QGate2.cpp
                     qgates/QControlledGate2.cpp
                     qgates/CRotationZ.cpp
                     qgates/CPhase.cpp
                     qgates/CX.cpp
                     qgates/CY.cpp
                     qgates/CZ.cpp
                     qgates/SWAP.cpp
                     qgates/iSWAP.cpp
                     io/QASMFile.cpp
                     io/util.cpp
           )
target_include_directories( qclabpp PUBLIC ${PROJECT_SOURCE_DIR}/include )
target_compile_features( qclabpp PUBLIC cxx_std_17 )

# openmp
if( TARGET OpenMP::OpenMP_CXX )
  target_link_libraries( qclabpp PUBLIC OpenMP::OpenMP_CXX )
endif()

# openmp offloading
option( QCLAB_OMP_OFFLOADING "Enable OMP offloading" OFF )
if ( QCLAB_OMP_OFFLOADING )
  target_compile_definitions( qclabpp PUBLIC QCLAB_OMP_OFFLOADING )
endif()

