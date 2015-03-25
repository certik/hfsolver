include(LibFindMacros)

libfind_include(libint2.h libint)
libfind_library(int2 libint)

set(LIBINT_LIBRARIES ${INT2_LIBRARY})
set(LIBINT_INCLUDE_DIRS ${LIBINT_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libint DEFAULT_MSG
    LIBINT_LIBRARY LIBINT_INCLUDE_DIR)

mark_as_advanced(LIBINT_INCLUDE_DIR INT2_LIBRARY)
