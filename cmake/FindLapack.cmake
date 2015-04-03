include(LibFindMacros)

libfind_library(lapack lapack)
libfind_library(blas lapack)
set(LAPACK_LIBRARIES ${LAPACK_LIBRARY} ${BLAS_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Lapack DEFAULT_MSG LAPACK_LIBRARIES)

mark_as_advanced(LAPACK_LIBRARY BLAS_LIBRARY)
