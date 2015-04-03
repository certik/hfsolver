include(LibFindMacros)

libfind_include(petsc.h petsc)
libfind_library(petsc petsc)

set(PETSC_LIBRARIES ${PETSC_LIBRARY})
set(PETSC_INCLUDE_DIRS ${PETSC_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Petsc DEFAULT_MSG
    PETSC_LIBRARIES PETSC_INCLUDE_DIRS)

mark_as_advanced(PETSC_INCLUDE_DIR PETSC_LIBRARY)
