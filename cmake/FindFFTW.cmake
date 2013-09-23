# Note:
# fftw3.h is a C interface (all precisions)
# fftw3q.f03 is a double precision Fortran interface
find_path(FFTW_INCLUDE_DIR fftw3.h $ENV{PYTHONHPC}/include NO_DEFAULT_PATH)
find_path(FFTW_INCLUDE_DIR fftw3.h)
find_library(FFTW_LIBRARY fftw3 $ENV{PYTHONHPC}/lib NO_DEFAULT_PATH)
find_library(FFTW_LIBRARY fftw3)

set(FFTW_LIBRARIES ${FFTW_LIBRARY})
set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})

if(EXISTS "${FFTW_INCLUDE_DIR}/fftw3.f03")
    set(FFTW_INCLUDE_DIR_HAS_FORTRAN "OK")
else()
    message("Include dir '${FFTW_INCLUDE_DIR}' does not contain fftw3.f03")
    set(FFTW_INCLUDE_DIR_HAS_FORTRAN)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG
    FFTW_LIBRARY FFTW_INCLUDE_DIR FFTW_INCLUDE_DIR_HAS_FORTRAN)

mark_as_advanced(FFTW_INCLUDE_DIR FFTW_LIBRARY)
