find_library(BLOPEX_LIB1 BLOPEX $ENV{PYTHONHPC}/lib NO_DEFAULT_PATH)
find_library(BLOPEX_LIB2 blopex_serial_double $ENV{PYTHONHPC}/lib NO_DEFAULT_PATH)
set(BLOPEX_LIBRARIES ${BLOPEX_LIB1} ${BLOPEX_LIB2})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BLOPEX DEFAULT_MSG BLOPEX_LIBRARIES)

mark_as_advanced(BLOPEX_LIB1 BLOPEX_LIB2)
