# Specify the version being used
CMAKE_MINIMUM_REQUIRED (VERSION 3.13)

# Name your project here and the language(s) used
PROJECT (MsFEM C CXX Fortran)
set(CMAKE_CXX_STANDARD 14)

#--- Find package with cmake script

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

OPTION(USE_MKL "Use MKL" ON)
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

FIND_PACKAGE(BLAS)
IF (BLAS_FOUND)
    SET (EXTRALIB ${EXTRALIB} ${BLAS_LIBRARIES})
    IF (USE_MKL)
        SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DWITH_MKL")
    ENDIF (USE_MKL)
ENDIF (BLAS_FOUND)

FIND_PACKAGE(LAPACK)
IF (LAPACK_FOUND)
   SET (EXTRALIB ${EXTRALIB} ${LAPACK_LIBRARIES})
ENDIF (LAPACK_FOUND)

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

#--- Include self-contained METIS-4.0.3
#
# An alternative way is
#FIND_LIBRARY(METIS_LIB NAMES metis libmetis)
find_package(METIS REQUIRED)
#
#FILE(GLOB METIS_SOURCES "metis-4.0.3/*.c")
#ADD_LIBRARY (metis4 STATIC ${METIS_SOURCES})
#SET_TARGET_PROPERTIES (metis4
#                PROPERTIES COMPILE_FLAGS "-I. -fPIC")
#
#--- Include self-contained ARPACK
FILE(GLOB ARPACK_SOURCES "ARPACK/*")
ADD_LIBRARY (myEigSolver STATIC ${ARPACK_SOURCES})
SET_TARGET_PROPERTIES (myEigSolver
                PROPERTIES COMPILE_FLAGS "-I. -fPIC")

#--- Include SparseSolver
FILE(GLOB SPARSESOLVER_SOURCES "SparseSolver.d/*")
ADD_LIBRARY (mySparseSolver STATIC ${SPARSESOLVER_SOURCES} FortranAdapter.h)
SET_TARGET_PROPERTIES (mySparseSolver
                PROPERTIES COMPILE_FLAGS "-I. -fPIC")

#--- Include math/finite elements directories
FILE(GLOB FECODE_SOURCES "Elements/*" "Math/*")
ADD_LIBRARY (myFECODE STATIC ${FECODE_SOURCES})
SET_TARGET_PROPERTIES (myFECODE
                PROPERTIES COMPILE_FLAGS "-I. -fPIC")

if (USE_MKL)
    INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}" 
    		"${MKL_INCLUDE_DIR}")
ELSE (USE_MKL)
    INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}")
ENDIF (USE_MKL)

#--- Specify the include directories
if (USE_MKL)
    INCLUDE_DIRECTORIES ("${METIS_INCLUDE_DIRS}"
                "${PROJECT_SOURCE_DIR}/ARPACK"
                "${PROJECT_SOURCE_DIR}/SparseSolver.d" 
                "${PROJECT_SOURCE_DIR}/Utilities.d"
                "${PROJECT_SOURCE_DIR}/Elements"
                "${PROJECT_SOURCE_DIR}/Math"
		"${MKL_INCLUDE_DIR}")
ELSE (USE_MKL)
    INCLUDE_DIRECTORIES ("${METIS_INCLUDE_DIRS}"
                "${PROJECT_SOURCE_DIR}/ARPACK"
                "${PROJECT_SOURCE_DIR}/SparseSolver.d" 
                "${PROJECT_SOURCE_DIR}/Utilities.d"
                "${PROJECT_SOURCE_DIR}/Elements"
                "${PROJECT_SOURCE_DIR}/Math"
		)
ENDIF (USE_MKL)

# Define executable and source files
ADD_EXECUTABLE(MsFEM test2D.cpp SparseSolver.d/FvCSpan.h)
TARGET_LINK_LIBRARIES(MsFEM mySparseSolver ${METIS_LIBRARY} myEigSolver myFECODE
              ${EXTRALIB})

