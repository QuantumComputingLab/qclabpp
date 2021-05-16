add_executable( run run.cpp
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
                    qgates/PointerGate1.cpp
                    qgates/QGate2.cpp
                    qgates/RotationXX.cpp
                    qgates/RotationYY.cpp
                    qgates/RotationZZ.cpp
                    qgates/QControlledGate2.cpp
                    qgates/CNOT.cpp
                    qgates/SWAP.cpp
                    qgates/CPhase.cpp
                    qgates/PointerGate2.cpp
              )
target_link_libraries( run PUBLIC qclab gtest )
target_include_directories( run PUBLIC ${CMAKE_SOURCE_DIR}/test )

add_executable( QAngle run.cpp QAngle.cpp )
target_link_libraries( QAngle PUBLIC qclab gtest )
target_include_directories( QAngle PUBLIC ${CMAKE_SOURCE_DIR}/test )

add_executable( QRotation run.cpp QRotation.cpp )
target_link_libraries( QRotation PUBLIC qclab gtest )
target_include_directories( QRotation PUBLIC ${CMAKE_SOURCE_DIR}/test )

add_executable( QCircuit run.cpp QCircuit.cpp )
target_link_libraries( QCircuit PUBLIC qclab gtest )
target_include_directories( QCircuit PUBLIC ${CMAKE_SOURCE_DIR}/test )

