add_library(imsi_lib)
target_sources(imsi_lib
        PUBLIC
        Element.h
        Mesh.h
        MeshUtils.h
        PRIVATE
        Mesh.cpp
        MeshUtils.cpp
)

target_include_directories(imsi_lib
        PUBLIC
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
        Kokkos::kokkoskernels Kokkos::kokkos
)

target_link_libraries(imsi_lib PUBLIC Kokkos::kokkoskernels Kokkos::kokkos)
