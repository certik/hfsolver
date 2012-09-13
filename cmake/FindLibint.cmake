find_path(LIBINT_INCLUDE_DIR libint2.h /usr/include/libint2 $ENV{SPKG_LOCAL}/include/libint2)
find_library(LIBINT_LIBRARY int2)

set(LIBINT_LIBRARIES ${LIBINT_LIBRARY} )
set(LIBINT_INCLUDE_DIRS ${LIBINT_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libint DEFAULT_MSG
    LIBINT_LIBRARY LIBINT_INCLUDE_DIR)

mark_as_advanced(LIBINT_INCLUDE_DIR LIBINT_LIBRARY)
