add_executable(
  Laplace1D.exe
  Laplace1D.cpp
)

target_include_directories(Laplace1D.exe
               PUBLIC
               tacho Kokkos::kokkoskernels Kokkos::kokkos
               )

target_link_libraries(Laplace1D.exe
  PUBLIC
  tacho Kokkos::kokkoskernels Kokkos::kokkos ${EXTRALIB}
)

#
# Example for testing SIMD acceleration
#
add_executable(
        simdPolar.exe
        simdPolar.cpp
)

target_include_directories(simdPolar.exe
        PUBLIC
        Kokkos::kokkos
)

target_link_libraries(simdPolar.exe
        PUBLIC
        Kokkos::kokkos ${EXTRALIB}
)

#
# Example for testing TeamPolicy
#
add_executable(
        TeamPolicy.exe
        team_policy_solution.cpp
)

target_include_directories(TeamPolicy.exe
        PUBLIC
        Kokkos::kokkos
)

target_link_libraries(TeamPolicy.exe
        PUBLIC
        Kokkos::kokkos ${EXTRALIB}
)

#
# CUDA assembly example (requires Kokkos with CUDA support)
#
if(Kokkos_ENABLE_CUDA)
  add_executable(
          cuda_assembly.exe
          cuda_assembly.cpp
  )

  target_include_directories(cuda_assembly.exe
          PUBLIC
          imsi_lib tacho Kokkos::kokkoskernels Kokkos::kokkos
  )

  target_link_libraries(cuda_assembly.exe
          PUBLIC
          imsi_lib tacho Kokkos::kokkoskernels Kokkos::kokkos ${EXTRALIB}
  )

  message(STATUS "CUDA assembly example enabled")
else()
  message(STATUS "CUDA assembly example disabled (Kokkos_ENABLE_CUDA not set)")
endif()

