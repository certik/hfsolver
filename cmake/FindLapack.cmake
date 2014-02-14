find_library(LAPACK_LIB lapack $ENV{PYTHONHPC}/lib NO_DEFAULT_PATH)
find_library(BLAS_LIB blas $ENV{PYTHONHPC}/lib NO_DEFAULT_PATH)
set(LAPACK_LIBRARIES ${LAPACK_LIB} ${BLAS_LIB})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Lapack DEFAULT_MSG LAPACK_LIBRARIES)

mark_as_advanced(LAPACK_LIB BLAS_LIB)
