find_library(MUMPS_LIB1 dmumps $ENV{PYTHONHPC}/lib NO_DEFAULT_PATH)
find_library(MUMPS_LIB2 mumps_common $ENV{PYTHONHPC}/lib NO_DEFAULT_PATH)
find_library(MUMPS_LIB3 pord $ENV{PYTHONHPC}/lib NO_DEFAULT_PATH)
set(MUMPS_LIBRARIES ${MUMPS_LIB1} ${MUMPS_LIB2} ${MUMPS_LIB3})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS DEFAULT_MSG MUMPS_LIBRARIES)

mark_as_advanced(MUMPS_LIB1 MUMPS_LIB2 MUMPS_LIB3)
