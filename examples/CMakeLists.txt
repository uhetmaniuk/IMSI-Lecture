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
        Kokkos::kokkos Vc
)

target_link_libraries(simdPolar.exe
        PUBLIC
        Kokkos::kokkos Vc ${EXTRALIB}
)
