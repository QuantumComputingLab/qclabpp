add_executable( qclab_tests run.cpp
                            QAngle.cpp
                            QRotation.cpp
                            QCircuit.cpp
                            dense/memory.cpp
                            dense/SquareMatrix.cpp
                            dense/transpose.cpp
                            dense/kron.cpp
                            qgates/QGate1.cpp
                            qgates/Identity.cpp
                            qgates/Hadamard.cpp
                            qgates/PauliX.cpp
                            qgates/PauliY.cpp
                            qgates/PauliZ.cpp
                            qgates/RotationX.cpp
                            qgates/RotationY.cpp
                            qgates/RotationZ.cpp
                            qgates/Phase.cpp
                            qgates/Phase45.cpp
                            qgates/Phase90.cpp
                            qgates/PointerGate1.cpp
                            qgates/QGate2.cpp
                            qgates/RotationXX.cpp
                            qgates/RotationYY.cpp
                            qgates/RotationZZ.cpp
                            qgates/QControlledGate2.cpp
                            qgates/CX.cpp
                            qgates/CY.cpp
                            qgates/CZ.cpp
                            qgates/CNOT.cpp
                            qgates/SWAP.cpp
                            qgates/iSWAP.cpp
                            qgates/CRotationX.cpp
                            qgates/CRotationY.cpp
                            qgates/CRotationZ.cpp
                            qgates/CPhase.cpp
                            qgates/PointerGate2.cpp
              )
target_link_libraries( qclab_tests PUBLIC qclabpp gtest )
target_include_directories( qclab_tests PUBLIC ${PROJECT_SOURCE_DIR}/test )

