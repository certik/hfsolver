set(SRC
    constants.f90   sorting.f90           types.f90
    utils.f90
    mesh.f90
    linalg.f90
    linalg_feast.f90
    lapack.f90
    feast.f90
    sparse.f90

    special_functions.f90
    solvers.f90
    hartree_screening.f90
    feutils.f90
    quadrature.f90
    quad_laguerre.f90
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
    fe_mesh.f90
    fourier.f90
    xc.f90

    splines.f90
    poisson3d_assembly.f90
    interp3d.f90
    isolve.f90

    c_hfsolver.f90

    optimize.f90

    random.f90

    integration.f90
    ofdft.f90
    ofdft_fft.f90
    ofdft_fe.f90
    md.f90
    bsplines.f90

    ewald.f90

    ffte.f90
)

# LIBINT
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

# PETSC
if(WITH_PETSC)
    set(SRC ${SRC}
        c_petsc.c
        petsc.f90
    )
else()
    set(SRC ${SRC}
        petsc.dum.f90
    )
endif()

# UMFPACK
if(WITH_UMFPACK)
    set(SRC ${SRC}
        umfpack.f90
    )
else()
    set(SRC ${SRC}
        umfpack.dum.f90
    )
endif()

# FFTW
if(WITH_FFTW)
    set(SRC ${SRC}
        fftw.f90
    )
endif()

# OPENMP
if(WITH_OPENMP)
    set(SRC ${SRC}
        openmp.f90
    )
else()
    set(SRC ${SRC}
        openmp.dum.f90
    )
endif()

if(WITH_MPI)
    set(SRC ${SRC}
        mpi2.f90
        mpi_dispatch.f90
        pffte.f90
        pofdft_fft.f90
    )
endif()

if(WITH_MKL)
    set(SRC ${SRC}
        c_mkl.c
        mkl.f90
    )
endif()

add_library(hfsolver ${SRC})

add_subdirectory(feast)
add_subdirectory(ffte)
add_subdirectory(tests)

target_link_libraries(hfsolver PUBLIC hfsolver::ffte)

install(TARGETS hfsolver
        RUNTIME DESTINATION bin
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
    )
