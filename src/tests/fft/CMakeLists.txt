macro(ADDPROG name)
    add_executable(${name} ${name}.f90)
    target_link_libraries(${name} hfsolver ${LIBS})
endmacro(ADDPROG)

project(fft)

addprog(sch1d_time_prop)
addtest(test_sch1d_osc)

if(WITH_MPI)
    ADDMPIEXECUTABLE(test_ffte_par)
    ADDMPITEST(test_ffte_par ffte1 1)
    ADDMPITEST(test_ffte_par ffte2 4)
    ADDMPITEST(test_ffte_par ffte3 4 2 2 2 32 32 32 1 2 2)
    ADDMPITEST(test_ffte_par ffte4 4 2 2 2 32 32 32 1 1 2)
    set_tests_properties(ffte4 PROPERTIES WILL_FAIL TRUE)
    ADDMPITEST(test_ffte_par ffte5 3 2 2 2 33 30 30 1 1 3)
    set_tests_properties(ffte5 PROPERTIES WILL_FAIL TRUE)
    ADDMPITEST(test_ffte_par ffte6 3 2 2 2 32 30 30 1 1 3)
    set_tests_properties(ffte6 PROPERTIES WILL_FAIL TRUE)
    ADDMPITEST(test_ffte_par ffte7 3 2 2 2 36 36 36 1 1 3)

    ADDMPIEXECUTABLE(test_ffte_par2)
    ADDMPITEST(test_ffte_par2 ffte_par2 4)
endif()

addtest(test_functional_derivative_min_fft)
set_tests_properties(test_functional_derivative_min_fft PROPERTIES LABELS slow)
