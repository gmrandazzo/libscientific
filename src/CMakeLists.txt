set(Lapack_SRCS
# Code for Matrix inversion and LU Decomposition
3dpart-lapack-3.8.0/dbdsdc.f
3dpart-lapack-3.8.0/dbdsqr.f
3dpart-lapack-3.8.0/dgebd2.f
3dpart-lapack-3.8.0/dgebrd.f
3dpart-lapack-3.8.0/dgelq2.f
3dpart-lapack-3.8.0/dgelqf.f
3dpart-lapack-3.8.0/dgeqr2.f
3dpart-lapack-3.8.0/dgeqrf.f
3dpart-lapack-3.8.0/dgesdd.f
3dpart-lapack-3.8.0/disnan.f
3dpart-lapack-3.8.0/dlabrd.f
3dpart-lapack-3.8.0/dlacpy.f
3dpart-lapack-3.8.0/dlaed6.f
3dpart-lapack-3.8.0/dlaisnan.f
3dpart-lapack-3.8.0/dlamch.f
3dpart-lapack-3.8.0/dlamrg.f
3dpart-lapack-3.8.0/dlange.f
3dpart-lapack-3.8.0/dlanst.f
3dpart-lapack-3.8.0/dlapy2.f
3dpart-lapack-3.8.0/dlarf.f
3dpart-lapack-3.8.0/dlarfb.f
3dpart-lapack-3.8.0/dlarfg.f
3dpart-lapack-3.8.0/dlarft.f
3dpart-lapack-3.8.0/dlartg.f
3dpart-lapack-3.8.0/dlas2.f
3dpart-lapack-3.8.0/dlascl.f
3dpart-lapack-3.8.0/dlasd0.f
3dpart-lapack-3.8.0/dlasd1.f
3dpart-lapack-3.8.0/dlasd2.f
3dpart-lapack-3.8.0/dlasd3.f
3dpart-lapack-3.8.0/dlasd4.f
3dpart-lapack-3.8.0/dlasd5.f
3dpart-lapack-3.8.0/dlasd6.f
3dpart-lapack-3.8.0/dlasd7.f
3dpart-lapack-3.8.0/dlasd8.f
3dpart-lapack-3.8.0/dlasda.f
3dpart-lapack-3.8.0/dlasdq.f
3dpart-lapack-3.8.0/dlasdt.f
3dpart-lapack-3.8.0/dlaset.f
3dpart-lapack-3.8.0/dlasq1.f
3dpart-lapack-3.8.0/dlasq2.f
3dpart-lapack-3.8.0/dlasq3.f
3dpart-lapack-3.8.0/dlasq4.f
3dpart-lapack-3.8.0/dlasq5.f
3dpart-lapack-3.8.0/dlasq6.f
3dpart-lapack-3.8.0/dlasr.f
3dpart-lapack-3.8.0/dlasrt.f
3dpart-lapack-3.8.0/dlassq.f
3dpart-lapack-3.8.0/dlasv2.f
3dpart-lapack-3.8.0/dorg2r.f
3dpart-lapack-3.8.0/dorgbr.f
3dpart-lapack-3.8.0/dorgl2.f
3dpart-lapack-3.8.0/dorglq.f
3dpart-lapack-3.8.0/dorgqr.f
3dpart-lapack-3.8.0/dorm2r.f
3dpart-lapack-3.8.0/dormbr.f
3dpart-lapack-3.8.0/dorml2.f
3dpart-lapack-3.8.0/dormlq.f
3dpart-lapack-3.8.0/dormqr.f
3dpart-lapack-3.8.0/ieeeck.f
3dpart-lapack-3.8.0/iladlc.f
3dpart-lapack-3.8.0/iladlr.f
3dpart-lapack-3.8.0/ilaenv.f
3dpart-lapack-3.8.0/iparmq.f
3dpart-lapack-3.8.0/lsame.f
3dpart-lapack-3.8.0/xerbla.f


3dpart-lapack-3.8.0/dgetri.f
3dpart-lapack-3.8.0/dgetrf.f
3dpart-lapack-3.8.0/dgetf2.f
3dpart-lapack-3.8.0/dlaswp.f
3dpart-lapack-3.8.0/dtrsm.f
3dpart-lapack-3.8.0/dtrtri.f
3dpart-lapack-3.8.0/dtrti2.f

3dpart-lapack-3.8.0/daxpy.f
3dpart-lapack-3.8.0/dbdsqr.f
3dpart-lapack-3.8.0/dcopy.f
3dpart-lapack-3.8.0/ddot.f
3dpart-lapack-3.8.0/dgebak.f
3dpart-lapack-3.8.0/dgebal.f
3dpart-lapack-3.8.0/dgebd2.f
3dpart-lapack-3.8.0/dgebrd.f
3dpart-lapack-3.8.0/dgeev.f
3dpart-lapack-3.8.0/dgehd2.f
3dpart-lapack-3.8.0/dgehrd.f
3dpart-lapack-3.8.0/dgelq2.f
3dpart-lapack-3.8.0/dgelqf.f
3dpart-lapack-3.8.0/dgemm.f
3dpart-lapack-3.8.0/dgemv.f
3dpart-lapack-3.8.0/dgeqr2.f
3dpart-lapack-3.8.0/dgeqrf.f
3dpart-lapack-3.8.0/dger.f
3dpart-lapack-3.8.0/dgesvd.f
3dpart-lapack-3.8.0/dhseqr.f
3dpart-lapack-3.8.0/disnan.f
3dpart-lapack-3.8.0/dlabad.f
3dpart-lapack-3.8.0/dlabrd.f
3dpart-lapack-3.8.0/dlacpy.f
3dpart-lapack-3.8.0/dladiv.f
3dpart-lapack-3.8.0/dlaexc.f
3dpart-lapack-3.8.0/dlahqr.f
3dpart-lapack-3.8.0/dlahr2.f
3dpart-lapack-3.8.0/dlaisnan.f
3dpart-lapack-3.8.0/dlaln2.f
3dpart-lapack-3.8.0/dlamch.f
3dpart-lapack-3.8.0/dlange.f
3dpart-lapack-3.8.0/dlanv2.f
3dpart-lapack-3.8.0/dlapy2.f
3dpart-lapack-3.8.0/dlaqr0.f
3dpart-lapack-3.8.0/dlaqr1.f
3dpart-lapack-3.8.0/dlaqr2.f
3dpart-lapack-3.8.0/dlaqr3.f
3dpart-lapack-3.8.0/dlaqr4.f
3dpart-lapack-3.8.0/dlaqr5.f
3dpart-lapack-3.8.0/dlarf.f
3dpart-lapack-3.8.0/dlarfb.f
3dpart-lapack-3.8.0/dlarfg.f
3dpart-lapack-3.8.0/dlarft.f
3dpart-lapack-3.8.0/dlarfx.f
3dpart-lapack-3.8.0/dlartg.f
3dpart-lapack-3.8.0/dlas2.f
3dpart-lapack-3.8.0/dlascl.f
3dpart-lapack-3.8.0/dlaset.f
3dpart-lapack-3.8.0/dlasq1.f
3dpart-lapack-3.8.0/dlasq2.f
3dpart-lapack-3.8.0/dlasq3.f
3dpart-lapack-3.8.0/dlasq4.f
3dpart-lapack-3.8.0/dlasq5.f
3dpart-lapack-3.8.0/dlasq6.f
3dpart-lapack-3.8.0/dlasr.f
3dpart-lapack-3.8.0/dlasrt.f
3dpart-lapack-3.8.0/dlassq.f
3dpart-lapack-3.8.0/dlasv2.f
3dpart-lapack-3.8.0/dlasy2.f
3dpart-lapack-3.8.0/dnrm2.f
3dpart-lapack-3.8.0/dorg2r.f
3dpart-lapack-3.8.0/dorgbr.f
3dpart-lapack-3.8.0/dorghr.f
3dpart-lapack-3.8.0/dorgl2.f
3dpart-lapack-3.8.0/dorglq.f
3dpart-lapack-3.8.0/dorgqr.f
3dpart-lapack-3.8.0/dorm2r.f
3dpart-lapack-3.8.0/dormbr.f
3dpart-lapack-3.8.0/dormhr.f
3dpart-lapack-3.8.0/dorml2.f
3dpart-lapack-3.8.0/dormlq.f
3dpart-lapack-3.8.0/dormqr.f
3dpart-lapack-3.8.0/drot.f
3dpart-lapack-3.8.0/dscal.f
3dpart-lapack-3.8.0/dswap.f
3dpart-lapack-3.8.0/dtrevc.f
3dpart-lapack-3.8.0/dtrexc.f
3dpart-lapack-3.8.0/dtrmm.f
3dpart-lapack-3.8.0/dtrmv.f
3dpart-lapack-3.8.0/idamax.f
3dpart-lapack-3.8.0/ieeeck.f
3dpart-lapack-3.8.0/iladlc.f
3dpart-lapack-3.8.0/iladlr.f
3dpart-lapack-3.8.0/ilaenv.f
3dpart-lapack-3.8.0/iparmq.f
3dpart-lapack-3.8.0/lsame.f
3dpart-lapack-3.8.0/lsamen.f
3dpart-lapack-3.8.0/xerbla.f
3dpart-lapack-3.8.0/dbdsqr.f
3dpart-lapack-3.8.0/dgetrf2.f
3dpart-lapack-3.8.0/dtrevc3.f)

set(Scientific_C_SRCS memwrapper.c
                   scientificinfo.c
                   optimization.c
                   vector.c
                   list.c
                   matrix.c
                   tensor.c
                   numeric.c
                   graphs.c
                   algebra.c
                   statistic.c
                   metricspace.c
                   interpolate.c
                   pca.c
                   pls.c
                   epls.c
                   #cpca.c
                   upca.c
                   upls.c
                   lda.c
#                    svd.c
                   clustering.c
                   mlr.c
                   modelvalidation.c)

set(Scientific_C_H
                scientificinfo.h
                optimization.h
                vector.h
                list.h
                matrix.h
                tensor.h
                numeric.h
                graphs.h
                algebra.h
                statistic.h
                metricspace.h
                interpolate.h
                pca.h
                pls.h
                epls.h
                #cpca.h
                upca.h
                upls.h
                lda.h
#                 svd.h
                clustering.h
                mlr.h
                modelvalidation.h)

set(NumericTest_SRC testnumeric.c)
set(StrVectorTest_SRC teststrvector.c)
set(DVectorTest_SRC testdvector.c)
set(UIVectorTest_SRC testuivector.c)
set(IVectorTest_SRC testivector.c)
set(VectorListTest_SRC testlist.c)
set(MatrixTest_SRC testmatrix.c)
set(TensorTest_SRC testtensor.c)
set(AlgebraTest_SRC testalgebra.c)
set(PCATest_SRC testpca.c)
set(PLSTest_SRC testpls.c)
set(EPLSTest_SRC testepls.c)
set(MPCATest_SRC testmwpca.c)
set(MPLSTest_SRC testmwpls.c)
set(CLUSTERINGTest_SRC testclustering.c)
set(MLRTest_SRC testmlr.c)
set(SVDTest_SRC testsvd.c)
set(CLASSIFIERTest_SRC testlda.c)
set(METRICSPACETest_SRC testmetrics.c)
set(INTERPOLATETest_SRC testinterpolate.c)
set(GRAPHSTest_SRC testgraphs.c)
set(OptimizationTest_SRC testoptimization.c)

if(WIN32)
  message(STATUS "Set mingw options")
  set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS_RELEASE)
  set(CMAKE_EXECUTABLE_SUFFIX ".exe")
  set(CMAKE_STATIC_LIBRARY_SUFFIX ".dll")
endif()

if(LAPACK_FOUND)
    add_library(libscientific SHARED ${Scientific_C_SRCS} ${Scientific_C_H})
    set_target_properties(libscientific PROPERTIES OUTPUT_NAME "scientific")
    target_link_libraries(libscientific m pthread ${LAPACK_LIBRARIES})
else()
    add_library(libscientific SHARED ${Scientific_C_SRCS} ${Scientific_C_H} ${Lapack_SRCS})
    set_target_properties(libscientific PROPERTIES OUTPUT_NAME "scientific")
    target_link_libraries(libscientific m pthread)
endif()

if(WIN32)
  install(TARGETS libscientific ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib COMPONENT library)
  install(FILES ${Scientific_C_H} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/scientific)
  install(FILES scientific.h DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
else()
  install(TARGETS libscientific LIBRARY DESTINATION ${CMAKE_INSTALL_FULL_LIBDIR})
  install(FILES ${Scientific_C_H} DESTINATION include/scientific)
  install(FILES scientific.h DESTINATION include)
endif()

add_executable(testnumeric ${NumericTest_SRC})
target_link_libraries(testnumeric libscientific)

add_executable(teststrvector ${StrVectorTest_SRC})
target_link_libraries(teststrvector libscientific)

add_executable(testdvector ${DVectorTest_SRC})
target_link_libraries(testdvector libscientific)

add_executable(testuivector ${UIVectorTest_SRC})
target_link_libraries(testuivector libscientific)

add_executable(testivector ${IVectorTest_SRC})
target_link_libraries(testivector libscientific)

add_executable(testlist ${VectorListTest_SRC})
target_link_libraries(testlist libscientific)

add_executable(testmatrix ${MatrixTest_SRC})
target_link_libraries(testmatrix libscientific)

add_executable(testtensor ${TensorTest_SRC})
target_link_libraries(testtensor libscientific)

add_executable(testalgebra ${AlgebraTest_SRC})
target_link_libraries(testalgebra libscientific)

add_executable(testpca ${PCATest_SRC})
target_link_libraries(testpca libscientific)

add_executable(testpls ${PLSTest_SRC})
target_link_libraries(testpls libscientific)

add_executable(testepls ${EPLSTest_SRC})
target_link_libraries(testepls libscientific)

add_executable(testmwpca ${MPCATest_SRC})
target_link_libraries(testmwpca libscientific)

add_executable(testmwpls ${MPLSTest_SRC})
target_link_libraries(testmwpls libscientific)

add_executable(testclustering ${CLUSTERINGTest_SRC})
target_link_libraries(testclustering libscientific)

add_executable(testmlr ${MLRTest_SRC})
target_link_libraries(testmlr libscientific)

add_executable(testlda ${CLASSIFIERTest_SRC})
target_link_libraries(testlda libscientific)

add_executable(testmetricspace ${METRICSPACETest_SRC})
target_link_libraries(testmetricspace libscientific)

add_executable(testinterpolate ${INTERPOLATETest_SRC})
target_link_libraries(testinterpolate libscientific)

add_executable(testgraphs ${GRAPHSTest_SRC})
target_link_libraries(testgraphs libscientific)

add_executable(testoptimization ${OptimizationTest_SRC})
target_link_libraries(testoptimization libscientific)
