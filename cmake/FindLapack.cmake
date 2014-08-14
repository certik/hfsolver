find_library(LIBLAPACK lapack "${LAPACK_DIR}/lib")
find_library(LIBBLAS blas "${LAPACK_DIR}/lib")
set(LAPACK_LIBRARIES ${LIBLAPACK} ${LIBBLAS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Lapack DEFAULT_MSG LAPACK_LIBRARIES)

mark_as_advanced(LIBLAPACK LIBBLAS)
