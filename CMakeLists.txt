# Specify the version being used
CMAKE_MINIMUM_REQUIRED (VERSION 3.23)

# Name your project here and the language(s) used
PROJECT(IMSI CXX Fortran)
set(CMAKE_CXX_STANDARD 17)

#--- Find package with cmake script

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

option(USE_MKL "Use MKL" OFF)
IF (USE_MKL)
   find_package(MKL REQUIRED)
   if (BLAS_mkl_core_LIBRARY)
       message("BLAS seems to be MKL")
       get_filename_component(MKL_LIB_Base ${BLAS_mkl_core_LIBRARY} DIRECTORY)
       get_filename_component(MKL_Base ${MKL_LIB_Base} DIRECTORY)
       find_path(MKL_INCLUDE_PATH "mkl.h" HINTS ${MKL_Base} PATH_SUFFIXES include)
       message("MKL_INCLUDE_PATH ${MKL_INCLUDE_PATH}")
   endif()
   SET(TPL_ENABLE_MKL ON)
ENDIF (USE_MKL)

option(USE_TACHO "Use Tacho sparse solver library" OFF)
if (USE_TACHO)
    SET(IMSI_USE_TACHO ON)
endif (USE_TACHO)

find_package(BLAS)
if (BLAS_FOUND)
    set (EXTRALIB ${EXTRALIB} ${BLAS_LIBRARIES})
    # Note: We use BLAS for Fortran sparse solver, but NOT for KokkosKernels
    # Do NOT set TPL_ENABLE_BLAS here as it would force KokkosKernels to enable BLAS
    if (USE_MKL)
        SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MKL")
    else()
        ENABLE_LANGUAGE(C)
        ENABLE_LANGUAGE(Fortran)
        INCLUDE(FortranCInterface)
        IF (FortranCInterface_GLOBAL_SUFFIX STREQUAL "")
          SET(F77_BLAS_MANGLE "(name,NAME) ${FortranCInterface_GLOBAL_PREFIX}name")
        ELSE ()
          SET(F77_BLAS_MANGLE "(name,NAME) ${FortranCInterface_GLOBAL_PREFIX}name ## ${FortranCInterface_GLOBAL_SUFFIX}")
        ENDIF ()
    endif(USE_MKL)
endif(BLAS_FOUND)

find_package(LAPACK)
IF (LAPACK_FOUND)
    SET (EXTRALIB ${EXTRALIB} ${LAPACK_LIBRARIES})
    # Note: We use LAPACK for Fortran sparse solver, but NOT for KokkosKernels
    # Do NOT set TPL_ENABLE_LAPACK here as it would force KokkosKernels to enable LAPACK
ENDIF (LAPACK_FOUND)

find_package(METIS)
IF (METIS_FOUND)
    SET (EXTRALIB ${EXTRALIB} ${METIS_LIBRARIES})
    # METIS is used by KokkosKernels for graph partitioning, so we set the TPL flag
    SET(TPL_ENABLE_METIS ON)
ENDIF (METIS_FOUND)

#--- KokkosKernels configuration
# Note: We use KokkosGraph (coloring) and KokkosSparse (transpose), but cannot disable
# BLAS/LAPACK/BATCHED components due to internal KokkosKernels dependencies.
# SPARSE requires BLAS headers even for basic graph operations.
message(STATUS "KokkosKernels: All components enabled (SPARSE/GRAPH require BLAS/LAPACK dependencies)")

add_subdirectory(TPL)
add_subdirectory(src)
add_subdirectory(examples)

###################

add_executable(main.exe  main.cpp)

if(USE_TACHO)
  target_include_directories(main.exe
          PUBLIC
          imsi_lib tacho Kokkos::kokkoskernels Kokkos::kokkos
  )

  target_link_libraries(main.exe
          PUBLIC
          imsi_lib tacho Kokkos::kokkoskernels Kokkos::kokkos ${EXTRALIB}
  )
else()
  target_include_directories(main.exe
          PUBLIC
          imsi_lib Kokkos::kokkoskernels Kokkos::kokkos
  )

  target_link_libraries(main.exe
          PUBLIC
          imsi_lib Kokkos::kokkoskernels Kokkos::kokkos ${EXTRALIB}
  )
endif()
