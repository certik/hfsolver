include(LibFindMacros)

libfind_library(dmumps mumps)
libfind_library(mumps_common mumps)
libfind_library(pord mumps)
set(MUMPS_LIBRARIES ${DMUMPS_LIBRARY} ${MUMPS_COMMON_LIBRARY} ${PORD_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS DEFAULT_MSG DMUMPS_LIBRARY
    MUMPS_COMMON_LIBRARY PORD_LIBRARY)

mark_as_advanced(DMUMPS_LIBRARY MUMPS_COMMON_LIBRARY PORD_LIBRARY)
