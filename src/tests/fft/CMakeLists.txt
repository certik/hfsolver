macro(ADDPROG name)
    add_executable(${name} ${name}.f90)
    target_link_libraries(${name} hfsolver ${LIBS})
endmacro(ADDPROG)

project(fft)

addprog(sch1d_time_prop)
addtest(test_sch1d_osc)
