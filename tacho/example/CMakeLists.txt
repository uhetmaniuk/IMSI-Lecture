
add_executable(
  Tacho_ExampleCombineDataFileToMatrixMarketFile.x
  Tacho_ExampleCombineDataFileToMatrixMarketFile.cpp
)

target_include_directories(Tacho_ExampleCombineDataFileToMatrixMarketFile.x
               PUBLIC
               tacho Kokkos::kokkos
               )

target_link_libraries(Tacho_ExampleCombineDataFileToMatrixMarketFile.x
             tacho)

#
# Driver double
#
add_executable(
  Tacho_ExampleDriverDouble.x
  Tacho_ExampleDriver_double.cpp
)

target_include_directories(Tacho_ExampleDriverDouble.x
               PUBLIC
               tacho Kokkos::kokkos
               )

target_link_libraries(Tacho_ExampleDriverDouble.x
             tacho)

#
# Driver float
#
add_executable(
  Tacho_ExampleDriverFloat.x
  Tacho_ExampleDriver_float.cpp
)

target_include_directories(Tacho_ExampleDriverFloat.x
               PUBLIC
               tacho Kokkos::kokkos
               )

target_link_libraries(Tacho_ExampleDriverFloat.x
             tacho)

#
# External Interface
# 
add_executable(
  Tacho_ExampleExternalInterface.x
  Tacho_ExampleExternalInterface.cpp
)

target_include_directories(Tacho_ExampleExternalInterface.x
               PUBLIC
               tacho Kokkos::kokkos
               )

target_link_libraries(Tacho_ExampleExternalInterface.x
             tacho)

#
# NVIDIA cuSolver
# 
IF(Kokkos_ENABLE_CUDA) 
  add_executable(
    Tacho_ExampleCholCuSolver.x
    Tacho_ExampleCholCuSolver.cpp
  )
  target_include_directories(Tacho_ExampleCholCuSolver.x
                 PUBLIC
                 tacho Kokkos::kokkos
                 )
  target_link_libraries(Tacho_ExampleCholCuSolver.x
               tacho)
ENDIF()

#
# Intel MKL Pardiso 
#
IF(Kokkos_ENABLE_OPENMP OR Kokkos_ENABLE_SERIAL)
  IF(USE_MKL)
    add_executable(
      Tacho_ExampleCholPardiso.x
      Tacho_ExampleCholPardiso.cpp
    )
    target_include_directories(Tacho_ExampleCholPardiso.x
                   PUBLIC
                   tacho Kokkos::kokkos
                   )
    target_link_libraries(Tacho_ExampleCholPardiso.x
                 tacho)
  ENDIF()
ENDIF()

