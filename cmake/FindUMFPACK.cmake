include(LibFindMacros)

libfind_include(umfpack.h umfpack)
libfind_library(umfpack umfpack)
libfind_library(cholmod umfpack)
libfind_library(camd umfpack)
libfind_library(ccolamd umfpack)
libfind_library(colamd umfpack)
libfind_library(amd umfpack)
libfind_library(suitesparseconfig umfpack)

set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARY} ${CHOLMOD_LIBRARY} ${CAMD_LIBRARY}
    ${CCOLAMD_LIBRARY} ${COLAMD_LIBRARY} ${AMD_LIBRARY}
    ${SUITESPARSECONFIG_LIBRARY})
set(UMFPACK_INCLUDE_DIRS ${UMFPACK_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(UMFPACK DEFAULT_MSG UMFPACK_LIBRARY
    CHOLMOD_LIBRARY CAMD_LIBRARY CCOLAMD_LIBRARY COLAMD_LIBRARY AMD_LIBRARY
    SUITESPARSECONFIG_LIBRARY UMFPACK_INCLUDE_DIR)

mark_as_advanced(UMFPACK_INCLUDE_DIR UMFPACK_LIBRARY CHOLMOD_LIBRARY
    CAMD_LIBRARY CCOLAMD_LIBRARY COLAMD_LIBRARY AMD_LIBRARY
    SUITESPARSECONFIG_LIBRARY)
