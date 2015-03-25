include(LibFindMacros)

libfind_library(bddcml bddcml)

set(BDDCML_LIBRARIES ${BDDCML_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BDDCML DEFAULT_MSG BDDCML_LIBRARIES)
mark_as_advanced(BDDCML_LIBRARY)
