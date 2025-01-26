add_library(imsi_lib)
target_sources(imsi_lib
        PUBLIC
        Element.h
        fe1DQ1.h
        fe1DQ2.h
        fe2DQ1.h
        fe2DQ2.h
        fe3DQ1.h
        fe3DQ2.h
        FunctionExamples.h
        Mesh.h
        MeshUtils.h
        QuadratureRule.h
        MathUtils.h
        Utils.h
        PRIVATE
        Mesh.cpp
        MeshUtils.cpp
        QuadratureRule.cpp
        Utils.cpp
)

target_include_directories(imsi_lib
        PUBLIC
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
        Kokkos::kokkoskernels Kokkos::kokkos
)

target_link_libraries(imsi_lib PUBLIC Kokkos::kokkoskernels Kokkos::kokkos)
