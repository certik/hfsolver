include(LibFindMacros)

libfind_library(BLOPEX blopex)
libfind_library(blopex_serial_double blopex)
set(BLOPEX_LIBRARIES ${BLOPEX_LIBRARY} ${BLOPEX_SERIAL_DOUBLE_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BLOPEX DEFAULT_MSG BLOPEX_LIBRARIES)

mark_as_advanced(BLOPEX_LIBRARY BLOPEX_SERIAL_DOUBLE_LIBRARY)
