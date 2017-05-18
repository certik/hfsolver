find_package(Lapack REQUIRED)


find_library(ARPACK_LIBRARY NAMES arpack)
find_library(PARPACK_LIBRARY NAMES parpack)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ARPACK DEFAULT_MSG ARPACK_LIBRARY
    PARPACK_LIBRARY)

add_library(p::arpack INTERFACE IMPORTED)
set_property(TARGET p::arpack PROPERTY INTERFACE_LINK_LIBRARIES
    ${ARPACK_LIBRARY} ${PARPACK_LIBRARY} ${LAPACK_LIBRARIES})
