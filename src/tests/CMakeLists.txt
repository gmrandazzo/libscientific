
include_directories (${PROJECT_SOURCE_DIR}/src)

set(TestScientificInfo_SRC testscientificinfo.c)
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
set(CPCATest_SRC testcpca.c)
set(PLSTest_SRC testpls.c)
set(ICATest_SRC testica.c)
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

add_executable(testscientificinfo ${TestScientificInfo_SRC})
target_link_libraries(testscientificinfo libscientific)
add_test(NAME testscientificinfo COMMAND testscientificinfo)

add_executable(testnumeric ${NumericTest_SRC})
target_link_libraries(testnumeric libscientific)
add_test(NAME testnumeric COMMAND testnumeric)

add_executable(teststrvector ${StrVectorTest_SRC})
target_link_libraries(teststrvector libscientific)
add_test(NAME teststrvector COMMAND teststrvector)

add_executable(testdvector ${DVectorTest_SRC})
target_link_libraries(testdvector libscientific)
add_test(NAME testdvector COMMAND testdvector)

add_executable(testuivector ${UIVectorTest_SRC})
target_link_libraries(testuivector libscientific)
add_test(NAME testuivector COMMAND testuivector)

add_executable(testivector ${IVectorTest_SRC})
target_link_libraries(testivector libscientific)
add_test(NAME testivector COMMAND testivector)

add_executable(testlist ${VectorListTest_SRC})
target_link_libraries(testlist libscientific)
add_test(NAME testlist COMMAND testlist)

add_executable(testmatrix ${MatrixTest_SRC})
target_link_libraries(testmatrix libscientific)
add_test(NAME testmatrix COMMAND testmatrix)

add_executable(testtensor ${TensorTest_SRC})
target_link_libraries(testtensor libscientific)
add_test(NAME testtensor COMMAND testtensor)

add_executable(testalgebra ${AlgebraTest_SRC})
target_link_libraries(testalgebra libscientific)
add_test(NAME testalgebra COMMAND testalgebra)

add_executable(testpca ${PCATest_SRC})
target_link_libraries(testpca libscientific)
add_test(NAME testpca COMMAND testpca)

add_executable(testcpca ${CPCATest_SRC})
target_link_libraries(testcpca libscientific)
add_test(NAME testcpca COMMAND testcpca)

add_executable(testpls ${PLSTest_SRC})
target_link_libraries(testpls libscientific)
add_test(NAME testpls COMMAND testpls)

add_executable(testepls ${EPLSTest_SRC})
target_link_libraries(testepls libscientific)
add_test(NAME testepls COMMAND testepls)

add_executable(testmwpca ${MPCATest_SRC})
target_link_libraries(testmwpca libscientific)
add_test(NAME testmwpca COMMAND testmwpca)

add_executable(testmwpls ${MPLSTest_SRC})
target_link_libraries(testmwpls libscientific)
add_test(NAME testmwpls COMMAND testmwpls)

add_executable(testclustering ${CLUSTERINGTest_SRC})
target_link_libraries(testclustering libscientific)
add_test(NAME testclustering COMMAND testclustering)

add_executable(testmlr ${MLRTest_SRC})
target_link_libraries(testmlr libscientific)
add_test(NAME testmlr COMMAND testmlr)

add_executable(testlda ${CLASSIFIERTest_SRC})
target_link_libraries(testlda libscientific)
add_test(NAME testlda COMMAND testlda)

add_executable(testmetricspace ${METRICSPACETest_SRC})
target_link_libraries(testmetricspace libscientific)
add_test(NAME testmetricspace COMMAND testmetricspace)

add_executable(testinterpolate ${INTERPOLATETest_SRC})
target_link_libraries(testinterpolate libscientific)
add_test(NAME testinterpolate COMMAND testinterpolate)

add_executable(testgraphs ${GRAPHSTest_SRC})
target_link_libraries(testgraphs libscientific)
add_test(NAME testgraphs COMMAND testgraphs)

add_executable(testoptimization ${OptimizationTest_SRC})
target_link_libraries(testoptimization libscientific)
add_test(NAME testoptimization COMMAND testoptimization)
