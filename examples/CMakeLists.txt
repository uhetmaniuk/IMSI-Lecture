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
  tacho Kokkos::kokkoskernels Kokkos::kokkos
)

#add_executable(
#        Laplace2D9pt.exe
#        Laplace2D9pt.cpp
#)
#
#target_include_directories(Laplace2D9pt.exe
#        PUBLIC
#        tacho Kokkos::kokkos
#)
#
#target_link_libraries(Laplace2D9pt.exe
#        PUBLIC
#        tacho
#)
