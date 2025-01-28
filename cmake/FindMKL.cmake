################################################################################
#
# \file      cmake/FindMKL.cmake
# \author    J. Bakosi
# \copyright 2012-2015, Jozsef Bakosi, 2016, Los Alamos National Security, LLC.
# \brief     Find the Math Kernel Library from Intel
# \date      Thu 26 Jan 2017 02:05:50 PM MST
#
################################################################################

# Find the Math Kernel Library from Intel
#
#  MKL_FOUND - System has MKL
#  MKL_INCLUDE_DIRS - MKL include files directories
#  MKL_LIBRARIES - The MKL libraries
#  MKL_INTERFACE_LIBRARY - MKL interface library
#  MKL_SEQUENTIAL_LAYER_LIBRARY - MKL sequential layer library
#  MKL_CORE_LIBRARY - MKL core library
#
#  The environment variables MKL_BASE and INTEL are used to find the library.
#  Everything else is ignored. If MKL is found "-DMKL_ILP64" is added to
#  CMAKE_C_FLAGS and CMAKE_CXX_FLAGS.
#
#  Example usage:
#
#  find_package(MKL)
#  if(MKL_FOUND)
#    target_link_libraries(TARGET ${MKL_LIBRARIES})
#  endif()
if(WIN32)
    if(CYGWIN)
        set(INTEL_ROOT "/cygdrive/c/Program\ Files\ \(x86\)/IntelSWTools/compilers_and_libraries/windows/"
                CACHE PATH "Folder contains Intel libs")
        set(INTEL_TL_ROOT "${INTEL_ROOT}")
        set(MKL_BASE "${INTEL_ROOT}/mkl" CACHE PATH "Folder at the root of Intel MKL")
    else(CYGWIN)
        #set(INTEL_ROOT "/opt/intel" CACHE PATH "Folder contains Intel libs")
        set(INTEL_ROOT "C:\\Program Files (x86)\\IntelSWTools\\compilers_and_libraries\\windows"
                CACHE PATH "Folder contains Intel libs")
        set(INTEL_TL_ROOT "${INTEL_ROOT}")
        set(MKL_BASE "${INTEL_ROOT}/mkl" CACHE PATH "Folder at the root of Intel MKL")
    endif(CYGWIN)
else(WIN32)
    set(INTEL_ROOT "/opt/intel" CACHE PATH "Folder contains Intel libs")
    set(INTEL_TL_ROOT "${INTEL_ROOT}")
    set(MKL_BASE "${INTEL_ROOT}/mkl" CACHE PATH "Folder at the root of Intel MKL")
endif(WIN32)

# If already in cache, be silent
#if (MKL_INCLUDE_DIRS AND MKL_LIBRARIES AND MKL_INTERFACE_LIBRARY AND
#        MKL_SEQUENTIAL_LAYER_LIBRARY AND MKL_CORE_LIBRARY)
#    set (MKL_FIND_QUIETLY TRUE)
#endif()

#if(NOT BUILD_SHARED_LIBS)
#	set(INT_LIB "libmkl_intel_ilp64.a")
#	set(SEQ_LIB "libmkl_sequential.a")
#	set(THR_LIB "libmkl_intel_thread.a")
#	set(COR_LIB "libmkl_core.a")
#else()
if(MKL_64)
    set(INT_LIB "mkl_intel_ilp64")
    find_library(MKL_INTERFACE_LIBRARY_I64
            NAMES "mkl_intel_ilp64"
            PATHS ${MKL_BASE}/lib
            ${MKL_BASE}/lib/intel64
            $ENV{INTEL}/mkl/lib/intel64
            NO_DEFAULT_PATH)
    set(MKL_INTERFACE_LIBRARY ${MKL_INTERFACE_LIBRARY_I64})
    set(MKL_64_FLAGS "-DMKL_ILP64")
#    message("64 bit MKL")
else(MKL_64)
    set(INT_LIB "mkl_intel_lp64")
    find_library(MKL_INTERFACE_LIBRARY_I32
            NAMES "mkl_intel_lp64"
            PATHS ${MKL_BASE}/lib
            ${MKL_BASE}/lib/intel64
            $ENV{INTEL}/mkl/lib/intel64
            NO_DEFAULT_PATH)
    set(MKL_INTERFACE_LIBRARY ${MKL_INTERFACE_LIBRARY_I32})
#    message ("32 bit MKL")
endif(MKL_64)
set(SEQ_LIB "mkl_sequential")
set(THR_LIB "mkl_intel_thread")
set(COR_LIB "mkl_core")
#endif()

set(MKL_BASE $ENV{MKL_BASE} CACHE PATH "Path the the root directory of Intel MKL library")

find_path(MKL_INCLUDE_DIR NAMES mkl.h HINTS ${MKL_BASE}/include)
#message("MKL include dir is ${MKL_INCLUDE_DIR}")

#message("MKL will look for ${INT_LIB}")


find_library(MKL_SEQUENTIAL_LAYER_LIBRARY
        NAMES ${SEQ_LIB}
        PATHS ${MKL_BASE}/lib
        ${MKL_BASE}/lib/intel64
        $ENV{INTEL}/mkl/lib/intel64
        NO_DEFAULT_PATH)
#message("Sequential: ${MKL_SEQUENTIAL_LAYER_LIBRARY} for ${SEQ_LIB} in ${MKL_BASE}/lib")

find_library(MKL_THREADED_LAYER_LIBRARY
        NAMES ${THR_LIB}
        PATHS ${MKL_BASE}/lib
        ${MKL_BASE}/lib/intel64
        $ENV{INTEL}/mkl/lib/intel64
        NO_DEFAULT_PATH)
#message("Threaded: ${MKL_THREADED_LAYER_LIBRARY} for ${SEQ_LIB} in ${MKL_BASE}/lib")

find_library(MKL_CORE_LIBRARY
        NAMES ${COR_LIB}
        PATHS ${MKL_BASE}/lib
        ${MKL_BASE}/lib/intel64
        $ENV{INTEL}/mkl/lib/intel64
        NO_DEFAULT_PATH)

set(MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIR})

############################ RTL layer ##########################
if(WIN32)
    set(MKL_RTL_LIBNAME libiomp5md)
else()
    set(MKL_RTL_LIBNAME iomp5)
endif()
find_library(MKL_RTL_LIBRARY ${MKL_RTL_LIBNAME}
        PATHS ${INTEL_RTL_ROOT}/lib
        /opt/intel//compilers_and_libraries_2018.2.164/mac/compiler/lib/
        /opt/intel//compilers_and_libraries/mac/lib/
        )
#message("RTL ${MKL_RTL_LIBRARY}")
#message("MKL Interface library ${MKL_INTERFACE_LIBRARY}")
if(MKL_THREAD)
#    message("MKL Is threaded.")
#        message ("MKL_RTL: ${MKL_RTL_LIBRARY}")
    set(MKL_LIBRARIES ${MKL_INTERFACE_LIBRARY} ${MKL_THREADED_LAYER_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_RTL_LIBRARY})
else(MKL_THREAD)
#    message("MKL Is not threaded.")
    set(MKL_LIBRARIES ${MKL_INTERFACE_LIBRARY} ${MKL_SEQUENTIAL_LAYER_LIBRARY} ${MKL_CORE_LIBRARY})
endif(MKL_THREAD)


if (MKL_INCLUDE_DIR AND
        MKL_INTERFACE_LIBRARY AND
        MKL_SEQUENTIAL_LAYER_LIBRARY AND
        MKL_CORE_LIBRARY)

    if (NOT DEFINED ENV{CRAY_PRGENVPGI} AND
            NOT DEFINED ENV{CRAY_PRGENVGNU} AND
            NOT DEFINED ENV{CRAY_PRGENVCRAY} AND
            NOT DEFINED ENV{CRAY_PRGENVINTEL})
        set(ABI "-m64")
    endif()

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MKL_64_FLAGS} ${ABI}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MKL_64_FLAGS} ${ABI}")

else()

    set(MKL_INCLUDE_DIRS "")
    set(MKL_LIBRARIES "")
    set(MKL_INTERFACE_LIBRARY "")
    set(MKL_SEQUENTIAL_LAYER_LIBRARY "")
    set(MKL_CORE_LIBRARY "")

endif()

#message("MKL_LIBRARIES  ${MKL_LIBRARIES} Flags ${MKL_64_FLAGS}")


# Handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
if(MKL_THREAD)
#    message("Using threaded MKL")
    FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_DIRS MKL_INTERFACE_LIBRARY MKL_THREADED_LAYER_LIBRARY MKL_CORE_LIBRARY)

    MARK_AS_ADVANCED(MKL_INCLUDE_DIRS MKL_LIBRARIES MKL_INTERFACE_LIBRARY MKL_THREADED_LAYER_LIBRARY MKL_CORE_LIBRARY)
else(MKL_THREAD)
#    message("Using sequential MKL")
  FIND_PACKAGE_HANDLE_STANDARD_ARGS(MKL DEFAULT_MSG MKL_LIBRARIES MKL_INCLUDE_DIRS MKL_INTERFACE_LIBRARY MKL_SEQUENTIAL_LAYER_LIBRARY MKL_CORE_LIBRARY)

  MARK_AS_ADVANCED(MKL_INCLUDE_DIRS MKL_LIBRARIES MKL_INTERFACE_LIBRARY MKL_SEQUENTIAL_LAYER_LIBRARY MKL_CORE_LIBRARY)
endif(MKL_THREAD)