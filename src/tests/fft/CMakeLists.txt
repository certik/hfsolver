macro(ADDPROG name)
    add_executable(${name} ${name}.f90)
    target_link_libraries(${name} hfsolver ${LIBS})
endmacro(ADDPROG)

project(fft)

addprog(sch1d_time_prop)
addprog(ofmd_rt)
addtest(test_sch1d_osc)
addtest(test_sch1d_osc_eigen)
addtest(test_sch3d_osc)
addtest(test_sch3d_osc_eigen)
addtest(test_sch3d_rt)
