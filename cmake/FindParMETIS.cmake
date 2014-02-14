find_library(PARMETIS_LIB1 parmetis $ENV{PYTHONHPC}/lib NO_DEFAULT_PATH)
find_library(PARMETIS_LIB2 metis $ENV{PYTHONHPC}/lib NO_DEFAULT_PATH)
set(PARMETIS_LIBRARIES ${PARMETIS_LIB1} ${PARMETIS_LIB2})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ParMETIS DEFAULT_MSG PARMETIS_LIBRARIES)

mark_as_advanced(PARMETIS_LIB1 PARMETIS_LIB2)
