macro(ADDPROG name)
    add_executable(${name} ${name}.f90)
    target_link_libraries(${name} hfsolver ${LIBS})
endmacro(ADDPROG)

project(fem)

add_library(converged_energies converged_energies.f90)
target_link_libraries(converged_energies hfsolver)

addprog(poisson2d)
addprog(laplace2d_eig)
addprog(schroed2d)
addprog(schroed3d)
addprog(conv_poisson3d)
addprog(ofmd)
addprog(ofmd_time_propagation)

if(WITH_BDDC)
    addprog(test_free_energy_bddc)
endif()

addtest(test_poisson3d)
addtest(test_free_energy converged_energies)
addtest(test_free_energy2 converged_energies)
addtest(test_free_energy3 converged_energies)
addtest(test_free_energy4 converged_energies)
set_tests_properties(test_free_energy4 PROPERTIES LABELS MD)
addtest(test_minimization)
addtest(test_fft_ofmd)
set_tests_properties(test_fft_ofmd PROPERTIES LABELS MD)
addtest(test_free_energy_fft converged_energies)
set_tests_properties(test_free_energy_fft PROPERTIES LABELS MD)
addtest(test_free_energy_fft2 converged_energies)
set_tests_properties(test_free_energy_fft2 PROPERTIES LABELS MD)

addtest(test_functional_derivative converged_energies)

set_tests_properties(test_free_energy4 PROPERTIES LABELS slow)
set_tests_properties(test_functional_derivative PROPERTIES LABELS slow)

addtest(test_functional_derivative_fft converged_energies)
set_tests_properties(test_functional_derivative_fft PROPERTIES LABELS slow)
