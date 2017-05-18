project(bsplines)

macro(ADDTEST2 name)
    add_executable(${name} ${name}.f90)
    target_link_libraries(${name} ${ARGN} hfsolver ${LIBS})
    add_test(${name} ${PROJECT_BINARY_DIR}/${name})
endmacro(ADDTEST2)

add_library(schroed_util schroed_util.f90)

addtest(test_bsplines)
addtest2(test_ppum schroed_util)
addtest(test_dirac)
addtest2(test_schroed schroed_util)
