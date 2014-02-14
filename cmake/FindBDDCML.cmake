find_library(BDDCML_LIBRARIES bddcml $ENV{PYTHONHPC}/lib NO_DEFAULT_PATH)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BDDCML DEFAULT_MSG BDDCML_LIBRARIES)
