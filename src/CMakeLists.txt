add_library(imsi_lib)
target_sources(imsi_lib
        PUBLIC
        debug.h
        Element.h
        fe1DQ1.h
        fe1DQ2.h
        fe2DQ1.h
        fe2DQ2.h
        fe3DQ1.h
        fe3DQ2.h
        FunctionExamples.h
        linkfc.h
        MathUtils.h
        Mesh.h
        MeshUtils.h
        my_sparse_solver.h
        QuadratureRule.h
        SparseMatrix.hpp
        stat.h
        SymmetricSparse.hpp
        Utils.h
        PRIVATE
        Mesh.cpp
        MeshUtils.cpp
        QuadratureRule.cpp
        Utils.cpp
        assmb.f		igathr.f	sassmb.f	zblkslvn.f
        betree.f	inpnv.f		sblkl2.f	symfc2.f	zblkslvp.f
        bfinit.f	epost2.f	invinv.f	sblkldl.f	symfct.f	zero.f
        blkl2.f		etordr.f	ldindx.f	sblkns.f	szero.f		zgecp.f
        blkldl.f	etpost.f	mmdelm.f	sblkslv.f	zassmb.f	zgers.f
        blkns.f		etree.f		mmdint.f	sblkslvn.f	zblkl2.f	zmmpyi.f
        blkslv.f	fcnthn.f	mmdnum.f	sfinit.f	zblkldl.f
        blkslvn.f	fntsiz.f	mmdupd.f	sgecp.f		zblkns.f	zzeromat.f
        btree2.f	mmpyi.f		sgers.f		zblkslv.f   sblkl2.f
        chordr.f	fsup1.f		slacol.f	zblkslv2.f
        dpotf2_ldl.f	dpotrf_ldl.f	spotf2_ldl.f	spotrf_ldl.f	zpotf2_ldl.f
        dgecp.f		fsup2.f		zblkslv3.f
        dgers.f		genmmd.f	ordmmd.f	zblkslv4.f
)

target_include_directories(imsi_lib
        PUBLIC
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}
        tacho Kokkos::kokkoskernels Kokkos::kokkos
)

target_link_libraries(imsi_lib PUBLIC tacho Kokkos::kokkoskernels Kokkos::kokkos)
