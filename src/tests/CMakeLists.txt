include_directories(${PROJECT_BINARY_DIR}/src)

macro(ADDTEST name)
    add_executable(${name} ${name}.f90)
    target_link_libraries(${name} hfsolver ${ARGN} ${LIBS})
    add_test(${name} ${PROJECT_BINARY_DIR}/${name})
endmacro(ADDTEST)

macro(ADDMPIEXECUTABLE name)
    add_executable(${name} ${name}.f90)
    target_link_libraries(${name} hfsolver ${ARGN} ${LIBS})
endmacro(ADDMPIEXECUTABLE)

macro(ADDMPITEST name testname NCORES)
    add_test(${testname} mpirun -n ${NCORES} ${PROJECT_BINARY_DIR}/${name}
        ${ARGN})
    set_tests_properties(${testname} PROPERTIES
        LABELS mpi
        PROCESSORS ${NCORES})
endmacro(ADDMPITEST)

macro(ADDPROG name)
    add_executable(${name} ${name}.f90)
    target_link_libraries(${name} hfsolver ${LIBS})
endmacro(ADDPROG)

project(tests)

ADDTEST(test_basis)
ADDTEST(test_basis_sto)
ADDTEST(test_gaunt)
ADDTEST(test_ijkl)
ADDTEST(test_sto)
ADDTEST(test_sto_hm)
ADDTEST(test_sto_he)
ADDTEST(test_sto_he2)
ADDTEST(test_sto_he3)
ADDTEST(test_sto_he4)
ADDTEST(test_sto_he5)
ADDTEST(test_sto_be)
ADDTEST(test_sto_be3)
set_tests_properties(test_sto_be3 PROPERTIES LABELS slow)
ADDTEST(test_sto_be4)
ADDTEST(test_sto_be5)
ADDTEST(test_sto_be6)
set_tests_properties(test_sto_be6 PROPERTIES LABELS slow)
ADDTEST(test_sto_be8)
set_tests_properties(test_sto_be8 PROPERTIES LABELS slow)
ADDTEST(test_sto_ne)
ADDTEST(test_sto_ne4)
ADDTEST(test_sto_ne5)
ADDTEST(test_sto_ne6)
ADDTEST(test_sto_ne7)
set_tests_properties(test_sto_ne7 PROPERTIES LABELS slow)
ADDTEST(test_sto_mg)
ADDTEST(test_sto_mg6)
set_tests_properties(test_sto_mg6 PROPERTIES LABELS slow)
ADDTEST(test_sto_xe)
ADDTEST(test_sto_ra)
ADDTEST(test_fe_he)
ADDTEST(test_fe_ne)
ADDTEST(test_fe_ra)
ADDTEST(test_fe_ra2)
ADDTEST(test_gaussians_atoms)
ADDTEST(test_gaussians_molecules)

ADDTEST(test_fm)
ADDTEST(test_bessel)

ADDTEST(test_sparse)

ADDTEST(test_fourier)
ADDTEST(test_ffte)

ADDTEST(test_brent)

ADDTEST(test_random)

ADDTEST(test_linregress)

ADDTEST(test_md)
set_tests_properties(test_md PROPERTIES LABELS MD)
ADDTEST(analyze_md)
set_tests_properties(analyze_md PROPERTIES DEPENDS test_md)
set_tests_properties(analyze_md PROPERTIES LABELS MD)

ADDTEST(test_ewald)
ADDTEST(test_ewald2)
ADDTEST(test_ewald3)
ADDTEST(test_ewald4)
ADDTEST(test_ewald5)
ADDTEST(test_ewald6)
ADDPROG(ewald_conv4)
ADDPROG(ewald_convN)

ADDPROG(test_ofmd_fft_fe)
ADDPROG(ofmd_fft)

ADDTEST(test_poisson_fft)


if(WITH_PETSC)
    ADDTEST(test_petsc)
endif()

if(WITH_FFTW)
    ADDTEST(test_fftw)
    add_executable(test_fftw_c test_fftw_c.c)
    target_link_libraries(test_fftw_c ${LIBS})
    add_test(test_fftw_c ${PROJECT_BINARY_DIR}/test_fftw_c)
endif()

if(WITH_MKL)
    ADDTEST(test_mkl)
endif()

if(WITH_UMFPACK)
    ADDTEST(test_umfpack)
    ADDTEST(test_umfpack_factorize)
    add_executable(test_umfpack_c test_umfpack_c.c)
    target_link_libraries(test_umfpack_c ${LIBS})
    add_test(test_umfpack_c ${PROJECT_BINARY_DIR}/test_umfpack_c)
endif()

if(WITH_ARPACK)
    ADDTEST(test_arpack)
endif()

add_subdirectory(thesis)
add_subdirectory(fem)
add_subdirectory(feast)
add_subdirectory(fft)
add_subdirectory(bsplines)
