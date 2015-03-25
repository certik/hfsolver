include(LibFindMacros)
# Note:
# fftw3.h is a C interface (all precisions)
# fftw3q.f03 is a double precision Fortran interface
libfind_include(fftw3.h fftw)
libfind_library(fftw3 fftw)

set(FFTW_LIBRARIES ${FFTW3_LIBRARY})
set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})

if(EXISTS "${FFTW_INCLUDE_DIR}/fftw3.f03")
    set(FFTW_INCLUDE_DIR_HAS_FORTRAN "OK")
else()
    message("Include dir '${FFTW_INCLUDE_DIR}' does not contain fftw3.f03")
    set(FFTW_INCLUDE_DIR_HAS_FORTRAN)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG
    FFTW_LIBRARIES FFTW_INCLUDE_DIRS FFTW_INCLUDE_DIR_HAS_FORTRAN)

mark_as_advanced(FFTW_INCLUDE_DIR FFTW3_LIBRARY)
