set(SRC
    constants.f90   sorting.f90           types.f90
    utils.f90
    mesh.f90
    linalg.f90
    lapack.f90
    openmp.f90

    special_functions.f90
    solvers.f90
    hartree_screening.f90
    feutils.f90
    quadrature.f90
    fe.f90
    sto.f90
    radialscf.f90
    hfprint.f90
    drivers.f90
    scf.f90
    cints.c
    qc.c
    qc.f90
    basis.f90
    basis_aux.f90
    mbpt.f90
    mbpt2.f90
    mbpt3.f90
    mbpt4.f90
    gf.f90
    gf2.f90
    debye.f90
    debye_potential.f90
    debye_potential_series.f90
)

if(WITH_LIBINT)
    set(SRC ${SRC}
        qc_libint.cpp
        gaussians_libint.f90
    )
else()
    set(SRC ${SRC}
        gaussians_hfsolver.f90
    )
endif()

set(LAPACK_LIBS
    lapack
    blas
    )

include_directories($ENV{SPKG_LOCAL}/include/libint2)

add_subdirectory(erd)

add_library(hfsolver ${SRC})

add_subdirectory(tests)
