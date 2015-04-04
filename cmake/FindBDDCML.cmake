include(LibFindMacros)

libfind_include(bddcml_interface_fortran.mod bddcml)
libfind_library(bddcml bddcml)

set(BDDCML_INCLUDE_DIRS ${BDDCML_INCLUDE_DIR})
set(BDDCML_LIBRARIES ${BDDCML_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BDDCML DEFAULT_MSG BDDCML_LIBRARY BDDCML_INCLUDE_DIR)
mark_as_advanced(BDDCML_LIBRARY BDDCML_INCLUDE_DIR)
