cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

set(CMAKE_USER_MAKE_RULES_OVERRIDE ${CMAKE_SOURCE_DIR}/cmake/UserOverride.cmake)

enable_language(Fortran)

project(hfsolver)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# Make sure that CMAKE_BUILD_TYPE is either Debug or Release:
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug
        CACHE STRING "Build type (Debug, Release)" FORCE)
endif ()
if (NOT (CMAKE_BUILD_TYPE STREQUAL "Debug" OR
        CMAKE_BUILD_TYPE STREQUAL "Release"))
    message("${CMAKE_BUILD_TYPE}")
    message(FATAL_ERROR "CMAKE_BUILD_TYPE must be one of: Debug, Release (current value: '${CMAKE_BUILD_TYPE}')")
endif ()

set(WITH_OPENMP yes
    CACHE BOOL "Build with OpenMP")


if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    # gfortran
    # Enable this if you want to check for single/double corruption (and use
    # the Debug build):
    #set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fdefault-real-8")

    # For profiling:
    #set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg -fno-omit-frame-pointer")
    #set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -g -pg")
    # Use it like this:
    # $ ./test_sto_he3
    # $ gprof ./test_sto_he3 | gprof2dot.py -n 0 -e 0 | dot -Tpng -o output.png

    if (WITH_OPENMP)
        set(CMAKE_Fortran_FLAGS_DEBUG
            "${CMAKE_Fortran_FLAGS_DEBUG}   -fopenmp")
        set(CMAKE_Fortran_FLAGS_RELEASE
            "${CMAKE_Fortran_FLAGS_RELEASE} -fopenmp")
        set(CMAKE_EXE_LINKER_FLAGS_DEBUG
            "${CMAKE_EXE_LINKER_FLAGS_DEBUG}   -fopenmp")
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE
            "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -fopenmp")
    endif ()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (WITH_OPENMP)
        set(CMAKE_Fortran_FLAGS_DEBUG
            "${CMAKE_Fortran_FLAGS_DEBUG}   -openmp")
        set(CMAKE_Fortran_FLAGS_RELEASE
            "${CMAKE_Fortran_FLAGS_RELEASE} -openmp")
        set(CMAKE_EXE_LINKER_FLAGS_DEBUG
            "${CMAKE_EXE_LINKER_FLAGS_DEBUG}   -openmp")
        set(CMAKE_EXE_LINKER_FLAGS_RELEASE
            "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -openmp")
    endif ()
endif ()

set(WITH_PYTHON no
    CACHE BOOL "Build with Python wrappers")

if (WITH_PYTHON)
    find_package(Python REQUIRED)
    find_package(NumPy REQUIRED)
    find_package(Cython REQUIRED)
    add_subdirectory(hfsolver)
endif ()

# ------------------------------------------------------------
# hfsolver dependencies:

# FEAST is included in src/feast
set(LIBS ${LIBS} feast)

# Lapack
find_package(Lapack REQUIRED)
set(LIBS ${LIBS} ${LAPACK_LIBRARIES})

# libint2
set(WITH_LIBINT no
    CACHE BOOL "Build with libint (http://sourceforge.net/p/libint/) support")

if(WITH_LIBINT)
    find_package(Libint REQUIRED)
    include_directories(${LIBINT_INCLUDE_DIRS})
    set(LIBS ${LIBS} ${LIBINT_LIBRARIES})
endif()

# fftw
set(WITH_FFTW no
    CACHE BOOL "Build with FFTW support")

if(WITH_FFTW)
    find_package(FFTW REQUIRED)
    include_directories(${FFTW_INCLUDE_DIRS})
    set(LIBS ${LIBS} ${FFTW_LIBRARIES})
endif()


set(WITH_MPI no
    CACHE BOOL "Build with MPI support")

# PETSc
set(WITH_PETSC no
    CACHE BOOL "Build with PETSc support")

if(WITH_PETSC)
    set(WITH_MPI yes)

    # Setup PETSc:
    find_package(Petsc REQUIRED)
    include_directories(${PETSC_INCLUDE_DIRS})
    set(LIBS ${LIBS} ${PETSC_LIBRARIES})
endif()

# MPI
if (WITH_MPI)
    find_package(MPI REQUIRED)
    include_directories(${MPI_INCLUDE_PATH})
    set(LIBS ${LIBS} ${MPI_LIBRARIES})
    set(CMAKE_C_FLAGS_DEBUG
        "${CMAKE_C_FLAGS_DEBUG}   ${MPI_COMPILE_FLAGS}")
    set(CMAKE_C_FLAGS_RELEASE
        "${CMAKE_C_FLAGS_RELEASE} ${MPI_COMPILE_FLAGS}")
    set(CMAKE_CXX_FLAGS_DEBUG
        "${CMAKE_C_FLAGS_DEBUG}   ${MPI_COMPILE_FLAGS}")
    set(CMAKE_CXX_FLAGS_RELEASE
        "${CMAKE_C_FLAGS_RELEASE} ${MPI_COMPILE_FLAGS}")
    set(CMAKE_Fortran_FLAGS_DEBUG
        "${CMAKE_Fortran_FLAGS_DEBUG}   ${MPI_COMPILE_FLAGS}")
    set(CMAKE_Fortran_FLAGS_RELEASE
        "${CMAKE_Fortran_FLAGS_RELEASE} ${MPI_COMPILE_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS_DEBUG
        "${CMAKE_EXE_LINKER_FLAGS_DEBUG}   ${MPI_LINK_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE
        "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${MPI_LINK_FLAGS}")
endif()

# ------------------------------------------------------------

enable_testing()

add_subdirectory(src)

message("\n")
message("Configuration results")
message("---------------------")
message("Fortran compiler: ${CMAKE_Fortran_COMPILER}")
message("C compiler      : ${CMAKE_C_COMPILER}")
message("C++ compiler    : ${CMAKE_CXX_COMPILER}")
message("Build type: ${CMAKE_BUILD_TYPE}")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message("Fortran compiler flags: ${CMAKE_Fortran_FLAGS_DEBUG}")
    message("C compiler flags      : ${CMAKE_C_FLAGS_DEBUG}")
    message("C++ compiler flags    : ${CMAKE_CXX_FLAGS_DEBUG}")
else ()
    message("Fortran compiler flags: ${CMAKE_Fortran_FLAGS_RELEASE}")
    message("C compiler flags      : ${CMAKE_C_FLAGS_RELEASE}")
    message("C++ compiler flags    : ${CMAKE_CXX_FLAGS_RELEASE}")
endif ()
message("LIBS: ${LIBS}")
message("Installation prefix: ${CMAKE_INSTALL_PREFIX}")
message("LAPACK_LIBRARIES: ${LAPACK_LIBRARIES}")
message("WITH_LIBINT: ${WITH_LIBINT}")
if(WITH_LIBINT)
    message("LIBINT_INCLUDE_DIRS: ${LIBINT_INCLUDE_DIRS}")
    message("LIBINT_LIBRARIES: ${LIBINT_LIBRARIES}")
endif()
message("WITH_FFTW: ${WITH_FFTW}")
if(WITH_FFTW)
    message("FFTW_INCLUDE_DIRS: ${FFTW_INCLUDE_DIRS}")
    message("FFTW_LIBRARIES: ${FFTW_LIBRARIES}")
endif()
message("WITH_PETSC: ${WITH_PETSC}")
if(WITH_PETSC)
    message("PETSC_INCLUDE_DIRS: ${PETSC_INCLUDE_DIRS}")
    message("PETSC_LIBRARIES: ${PETSC_LIBRARIES}")
endif()
message("WITH_PYTHON: ${WITH_PYTHON}")
message("WITH_OPENMP: ${WITH_OPENMP}")
message("WITH_MPI: ${WITH_MPI}")
if (WITH_MPI)
    message("MPI_INCLUDE_PATH: ${MPI_INCLUDE_PATH}")
    message("MPI_LIBRARIES: ${MPI_LIBRARIES}")
    message("MPI_COMPILE_FLAGS: ${MPI_COMPILE_FLAGS}")
    message("MPI_LINK_FLAGS: ${MPI_LINK_FLAGS}")
endif()
