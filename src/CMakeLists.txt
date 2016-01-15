set(Lapack_SRCS
# Code for Matrix inversion and LU Decomposition
3dpart-lapack-3.5.0/dgetri.f
3dpart-lapack-3.5.0/dgetrf.f
3dpart-lapack-3.5.0/dgetf2.f
3dpart-lapack-3.5.0/dlaswp.f
3dpart-lapack-3.5.0/dtrsm.f
3dpart-lapack-3.5.0/dtrtri.f
3dpart-lapack-3.5.0/dtrti2.f

3dpart-lapack-3.5.0/daxpy.f
3dpart-lapack-3.5.0/dbdsqr.f
3dpart-lapack-3.5.0/dcopy.f
3dpart-lapack-3.5.0/ddot.f
3dpart-lapack-3.5.0/dgebak.f
3dpart-lapack-3.5.0/dgebal.f
3dpart-lapack-3.5.0/dgebd2.f
3dpart-lapack-3.5.0/dgebrd.f
3dpart-lapack-3.5.0/dgeev.f
3dpart-lapack-3.5.0/dgehd2.f
3dpart-lapack-3.5.0/dgehrd.f
3dpart-lapack-3.5.0/dgelq2.f
3dpart-lapack-3.5.0/dgelqf.f
3dpart-lapack-3.5.0/dgemm.f
3dpart-lapack-3.5.0/dgemv.f
3dpart-lapack-3.5.0/dgeqr2.f
3dpart-lapack-3.5.0/dgeqrf.f
3dpart-lapack-3.5.0/dger.f
3dpart-lapack-3.5.0/dgesvd.f
3dpart-lapack-3.5.0/dhseqr.f
3dpart-lapack-3.5.0/disnan.f
3dpart-lapack-3.5.0/dlabad.f
3dpart-lapack-3.5.0/dlabrd.f
3dpart-lapack-3.5.0/dlacpy.f
3dpart-lapack-3.5.0/dladiv.f
3dpart-lapack-3.5.0/dlaexc.f
3dpart-lapack-3.5.0/dlahqr.f
3dpart-lapack-3.5.0/dlahr2.f
3dpart-lapack-3.5.0/dlaisnan.f
3dpart-lapack-3.5.0/dlaln2.f
3dpart-lapack-3.5.0/dlamch.f
3dpart-lapack-3.5.0/dlange.f
3dpart-lapack-3.5.0/dlanv2.f
3dpart-lapack-3.5.0/dlapy2.f
3dpart-lapack-3.5.0/dlaqr0.f
3dpart-lapack-3.5.0/dlaqr1.f
3dpart-lapack-3.5.0/dlaqr2.f
3dpart-lapack-3.5.0/dlaqr3.f
3dpart-lapack-3.5.0/dlaqr4.f
3dpart-lapack-3.5.0/dlaqr5.f
3dpart-lapack-3.5.0/dlarf.f
3dpart-lapack-3.5.0/dlarfb.f
3dpart-lapack-3.5.0/dlarfg.f
3dpart-lapack-3.5.0/dlarft.f
3dpart-lapack-3.5.0/dlarfx.f
3dpart-lapack-3.5.0/dlartg.f
3dpart-lapack-3.5.0/dlas2.f
3dpart-lapack-3.5.0/dlascl.f
3dpart-lapack-3.5.0/dlaset.f
3dpart-lapack-3.5.0/dlasq1.f
3dpart-lapack-3.5.0/dlasq2.f
3dpart-lapack-3.5.0/dlasq3.f
3dpart-lapack-3.5.0/dlasq4.f
3dpart-lapack-3.5.0/dlasq5.f
3dpart-lapack-3.5.0/dlasq6.f
3dpart-lapack-3.5.0/dlasr.f
3dpart-lapack-3.5.0/dlasrt.f
3dpart-lapack-3.5.0/dlassq.f
3dpart-lapack-3.5.0/dlasv2.f
3dpart-lapack-3.5.0/dlasy2.f
3dpart-lapack-3.5.0/dnrm2.f
3dpart-lapack-3.5.0/dorg2r.f
3dpart-lapack-3.5.0/dorgbr.f
3dpart-lapack-3.5.0/dorghr.f
3dpart-lapack-3.5.0/dorgl2.f
3dpart-lapack-3.5.0/dorglq.f
3dpart-lapack-3.5.0/dorgqr.f
3dpart-lapack-3.5.0/dorm2r.f
3dpart-lapack-3.5.0/dormbr.f
3dpart-lapack-3.5.0/dormhr.f
3dpart-lapack-3.5.0/dorml2.f
3dpart-lapack-3.5.0/dormlq.f
3dpart-lapack-3.5.0/dormqr.f
3dpart-lapack-3.5.0/drot.f
3dpart-lapack-3.5.0/dscal.f
3dpart-lapack-3.5.0/dswap.f
3dpart-lapack-3.5.0/dtrevc.f
3dpart-lapack-3.5.0/dtrexc.f
3dpart-lapack-3.5.0/dtrmm.f
3dpart-lapack-3.5.0/dtrmv.f
3dpart-lapack-3.5.0/idamax.f
3dpart-lapack-3.5.0/ieeeck.f
3dpart-lapack-3.5.0/iladlc.f
3dpart-lapack-3.5.0/iladlr.f
3dpart-lapack-3.5.0/ilaenv.f
3dpart-lapack-3.5.0/iparmq.f
3dpart-lapack-3.5.0/lsame.f
3dpart-lapack-3.5.0/lsamen.f
3dpart-lapack-3.5.0/xerbla.f
3dpart-lapack-3.5.0/dbdsqr.f)

set(Scientific_C_SRCS memwrapper.c
                   scientificinfo.c
                   vector.c
                   matrix.c
                   array.c
                   numeric.c
                   graphs.c
                   algebra.c
                   statistic.c
                   metricspace.c
                   pca.c
                   pls.c
                   upca.c
                   upls.c
                   lda.c
#                    svd.c
                   clustering.c
                   mlr.c)

set(Scientific_C_H
                scientificinfo.h
                vector.h
                matrix.h
                array.h
                numeric.h
                graphs.h
                algebra.h
                statistic.h
                metricspace.h
                pca.h
                pls.h
                upca.h
                upls.h
                lda.h
#                 svd.h
                clustering.h
                mlr.h)

set(StrVectorTest_SRC teststrvector.c)
set(DVectorTest_SRC testdvector.c)
set(UIVectorTest_SRC testuivector.c)
set(IVectorTest_SRC testivector.c)
set(MatrixTest_SRC testmatrix.c)
set(ArrayTest_SRC testarray.c)
set(AlgebraTest_SRC testalgebra.c)
set(PCATest_SRC testpca.c)
set(PLSTest_SRC testpls.c)
set(MPCATest_SRC testmwpca.c)
set(MPLSTest_SRC testmwpls.c)
set(CLUSTERINGTest_SRC testclustering.c)
set(MLRTest_SRC testmlr.c)
set(SVDTest_SRC testsvd.c)
set(CLASSIFIERTest_SRC testlda.c)
set(METRICSPACETest_SRC testmetrics.c)
set(GRAPHSTest_SRC testgraphs.c)

if(WIN32)
  message(STATUS "Set mingw options")
  set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS_RELEASE)
  set(CMAKE_EXECUTABLE_SUFFIX ".exe")
  set(CMAKE_STATIC_LIBRARY_SUFFIX ".dll")
endif()

#find_library(STATICLIBGCC libgcc.a PATH  /usr/local/Cellar/gcc/4.9.1/lib/gcc/x86_64-apple-darwin14.0.0/4.9.1/)
#find_library(STATICLIBQUADMATH libquadmath.a PATH  /usr/local/Cellar/gcc/4.9.1/lib/gcc/x86_64-apple-darwin14.0.0/4.9.1/)
#find_library(STATICGFORTRAN libgfortran.a PATH /usr/local/Cellar/gcc/4.9.1/lib/gcc/x86_64-apple-darwin14.0.0/4.9.1/)

# IF EMSCRIPTEN COMPILE REMOVE SHARED!
add_library(libscientific SHARED ${Scientific_C_SRCS} ${Scientific_C_H} ${Lapack_SRCS})
set_target_properties(libscientific PROPERTIES OUTPUT_NAME "scientific")
#target_link_libraries(libscientific m pthread ${STATICLIBQUADMATH} ${STATICLIBGCC} ${STATICGFORTRAN})
target_link_libraries(libscientific m pthread)

install(TARGETS libscientific LIBRARY DESTINATION ${LIB_INSTALL_DIR}
                      ARCHIVE DESTINATION ${LIB_INSTALL_DIR}
                      RUNTIME DESTINATION bin)

install(FILES ${Scientific_C_H} DESTINATION ${INCLUDE_INSTALL_DIR}/scientific)
install(FILES scientific.h DESTINATION ${INCLUDE_INSTALL_DIR})

#if(cmake_build_type_tolower STREQUAL "debug")
  add_executable(teststrvector ${StrVectorTest_SRC})
  target_link_libraries(teststrvector libscientific)

  add_executable(testdvector ${DVectorTest_SRC})
  target_link_libraries(testdvector libscientific)

  add_executable(testuivector ${UIVectorTest_SRC})
  target_link_libraries(testuivector libscientific)

  add_executable(testivector ${IVectorTest_SRC})
  target_link_libraries(testivector libscientific)

  add_executable(testmatrix ${MatrixTest_SRC})
  target_link_libraries(testmatrix libscientific)

  add_executable(testarray ${ArrayTest_SRC})
  target_link_libraries(testarray libscientific)

  add_executable(testalgebra ${AlgebraTest_SRC})
  target_link_libraries(testalgebra libscientific)

  add_executable(testpca ${PCATest_SRC})
  target_link_libraries(testpca libscientific)

  add_executable(testpls ${PLSTest_SRC})
  target_link_libraries(testpls libscientific)

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

  add_executable(testgraphs ${GRAPHSTest_SRC})
  target_link_libraries(testgraphs libscientific)

#   add_executable(testsvd ${SVDTest_SRC})
#   target_link_libraries(testsvd libscientific)
#endif()
