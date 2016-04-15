macro(ADDPROG name)
    add_executable(${name} ${name}.f90)
    target_link_libraries(${name} hfsolver ${LIBS})
endmacro(ADDPROG)

addprog(sch1d_time_prop)
