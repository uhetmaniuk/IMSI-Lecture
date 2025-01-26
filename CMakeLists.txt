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
ENDIF (USE_MKL)

find_package(BLAS)
if (BLAS_FOUND)
    set (EXTRALIB ${EXTRALIB} ${BLAS_LIBRARIES})
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
ENDIF (LAPACK_FOUND)

add_subdirectory(TPL)
add_subdirectory(src)
add_subdirectory(examples)

###################

add_executable(main.exe  main.cpp)

target_include_directories(main.exe
        PUBLIC
        imsi_lib tacho Kokkos::kokkoskernels Kokkos::kokkos Eigen
)

target_link_libraries(main.exe
        PUBLIC
        imsi_lib tacho Kokkos::kokkoskernels Kokkos::kokkos
)
