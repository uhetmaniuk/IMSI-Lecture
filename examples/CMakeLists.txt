add_executable(
  Laplace1D.exe
  Laplace1D.cpp
)

target_include_directories(Laplace1D.exe
               PUBLIC
               tacho Kokkos::kokkos
               )

target_link_libraries(Laplace1D.exe
  PUBLIC
  tacho
)

