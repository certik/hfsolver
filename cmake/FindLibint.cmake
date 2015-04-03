include(LibFindMacros)

libfind_include(libint2.h libint)
if(LIBINT_INCLUDE_DIR STREQUAL "LIBINT_INCLUDE_DIR-NOTFOUND")
    # Older versions of the library had the libint2.h inside the libint2
    # subdirectory, so we try to find it there as well
    libfind_include(libint2/libint2.h libint)
endif()

libfind_library(int2 libint)

set(LIBINT_LIBRARIES ${INT2_LIBRARY})
# Libint requires to include both the the libint2 subdirectory and its parent
# directory
set(LIBINT_INCLUDE_DIRS ${LIBINT_INCLUDE_DIR} ${LIBINT_INCLUDE_DIR}/libint2)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libint DEFAULT_MSG
    LIBINT_LIBRARIES LIBINT_INCLUDE_DIR)

mark_as_advanced(LIBINT_INCLUDE_DIR INT2_LIBRARY)
