# Laplace1D example requires Tacho
if(USE_TACHO)
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

  message(STATUS "Laplace1D example enabled (requires Tacho)")
else()
  message(STATUS "Laplace1D example disabled (requires USE_TACHO=ON)")
endif()

#
# Example comparing different sparse linear solvers
#
add_executable(
        Solver_Comparison.exe
        Solver_Comparison.cpp
)

target_include_directories(Solver_Comparison.exe
        PUBLIC
        imsi_lib Kokkos::kokkos
)

target_link_libraries(Solver_Comparison.exe
        PUBLIC
        imsi_lib Kokkos::kokkos ${EXTRALIB}
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

  if(USE_TACHO)
    target_include_directories(cuda_assembly.exe
            PUBLIC
            imsi_lib tacho Kokkos::kokkoskernels Kokkos::kokkos
    )

    target_link_libraries(cuda_assembly.exe
            PUBLIC
            imsi_lib tacho Kokkos::kokkoskernels Kokkos::kokkos ${EXTRALIB}
    )
  else()
    target_include_directories(cuda_assembly.exe
            PUBLIC
            imsi_lib Kokkos::kokkoskernels Kokkos::kokkos
    )

    target_link_libraries(cuda_assembly.exe
            PUBLIC
            imsi_lib Kokkos::kokkoskernels Kokkos::kokkos ${EXTRALIB}
    )
  endif()

  message(STATUS "CUDA assembly example enabled")
else()
  message(STATUS "CUDA assembly example disabled (Kokkos_ENABLE_CUDA not set)")
endif()

#
# OpenMP assembly example (requires Kokkos with OpenMP support)
# Uses ScaledLaplacian with Kokkos::OpenMP execution space
#
if(Kokkos_ENABLE_OPENMP)
  add_executable(
          openmp_assembly.exe
          openmp_assembly.cpp
  )

  if(USE_TACHO)
    target_include_directories(openmp_assembly.exe
            PUBLIC
            imsi_lib tacho Kokkos::kokkoskernels Kokkos::kokkos
    )

    target_link_libraries(openmp_assembly.exe
            PUBLIC
            imsi_lib tacho Kokkos::kokkoskernels Kokkos::kokkos ${EXTRALIB}
    )
  else()
    target_include_directories(openmp_assembly.exe
            PUBLIC
            imsi_lib Kokkos::kokkoskernels Kokkos::kokkos
    )

    target_link_libraries(openmp_assembly.exe
            PUBLIC
            imsi_lib Kokkos::kokkoskernels Kokkos::kokkos ${EXTRALIB}
    )
  endif()

  message(STATUS "OpenMP assembly example enabled")
else()
  message(STATUS "OpenMP assembly example disabled (Kokkos_ENABLE_OPENMP not set)")
endif()

#
# MFEM assembly example with OpenMP (requires Kokkos with OpenMP support)
# Uses ScaledLaplacian with MFEM_L elements and Kokkos::OpenMP execution space
#
if(Kokkos_ENABLE_OPENMP)
  add_executable(
          mfem_assembly.exe
          mfem_assembly.cpp
  )

  if(USE_TACHO)
    target_include_directories(mfem_assembly.exe
            PUBLIC
            imsi_lib tacho Kokkos::kokkoskernels Kokkos::kokkos
    )

    target_link_libraries(mfem_assembly.exe
            PUBLIC
            imsi_lib tacho Kokkos::kokkoskernels Kokkos::kokkos ${EXTRALIB}
    )
  else()
    target_include_directories(mfem_assembly.exe
            PUBLIC
            imsi_lib Kokkos::kokkoskernels Kokkos::kokkos
    )

    target_link_libraries(mfem_assembly.exe
            PUBLIC
            imsi_lib Kokkos::kokkoskernels Kokkos::kokkos ${EXTRALIB}
    )
  endif()

  message(STATUS "MFEM assembly example enabled")
else()
  message(STATUS "MFEM assembly example disabled (Kokkos_ENABLE_OPENMP not set)")
endif()

