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

