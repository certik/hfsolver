include(LibFindMacros)

libfind_library(scalapack scalapack)
set(SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCALAPACK DEFAULT_MSG SCALAPACK_LIBRARIES)
mark_as_advanced(SCALAPACK_LIBRARY)
