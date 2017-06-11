project(bsplines)

macro(ADDTEST2 name)
    add_executable(${name} ${name}.f90)
    target_link_libraries(${name} ${ARGN} hfsolver ${LIBS})
    add_test(${name} ${PROJECT_BINARY_DIR}/${name})
endmacro(ADDTEST2)

add_library(radial_util radial_util.f90)
target_link_libraries(radial_util hfsolver ${LIBS} lapack)

addtest(test_bsplines)
addtest2(test_ppum radial_util)
addtest(test_dirac radial_util)
addtest(test_dirac_ppum radial_util)
addtest2(test_schroed radial_util)
