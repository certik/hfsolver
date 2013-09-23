module fftw
use iso_c_binding, only: c_int, c_intptr_t, c_ptr, c_int32_t, &
    c_double_complex, c_double, c_funptr, c_size_t, c_float_complex, c_float, &
    c_char
implicit none

! FFTW provides the correct interface, we just need to include it and create a
! proper Fortran module from it:
include 'fftw3.f03'  ! general definitions, double precision interface

end module
