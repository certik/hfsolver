find_path(PETSC_INCLUDE_DIR petsc.h /usr/include $ENV{SPKG_LOCAL}/include)
find_library(PETSC_LIBRARY petsc)

set(PETSC_LIBRARIES ${PETSC_LIBRARY} )
set(PETSC_INCLUDE_DIRS ${PETSC_INCLUDE_DIR} /usr/include/mpi)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Petsc DEFAULT_MSG
    PETSC_LIBRARY PETSC_INCLUDE_DIR)

mark_as_advanced(PETSC_INCLUDE_DIR PETSC_LIBRARY)
