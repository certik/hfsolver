macro(ADDPROG name)
    add_executable(${name} ${name}.f90)
    target_link_libraries(${name} hfsolver ${LIBS})
endmacro(ADDPROG)

project(fem)

addprog(poisson2d)
addprog(laplace2d_eig)
addprog(schroed2d)
addprog(schroed3d)
addprog(conv_poisson3d)
addprog(ofmd)

addtest(test_poisson3d)
addtest(test_free_energy)
