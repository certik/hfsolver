include(LibFindMacros)

libfind_library(parmetis parmetis)
libfind_library(metis parmetis)
set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ParMETIS DEFAULT_MSG PARMETIS_LIBRARIES)

mark_as_advanced(PARMETIS_LIBRARIES METIS_LIBRARY)
