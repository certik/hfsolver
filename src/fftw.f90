module fftw
! Needed for fftw.f03:
use iso_c_binding, only: c_int, c_intptr_t, c_ptr, c_int32_t, &
    c_double_complex, c_double, c_funptr, c_size_t, c_float_complex, c_float, &
    c_char, c_loc, c_f_pointer
! Needed for fftw3l.f03:
use iso_c_binding, only: c_long_double_complex, c_long_double
use types, only: dp
implicit none

! FFTW provides the correct interface, we just need to include it and create a
! proper Fortran module from it:
include 'fftw3.f03'  ! general definitions, double precision interface

contains

! Helper functions for allocating/deallocating FFTW aligned arrays

function alloc1d(n1) result(x)
integer, intent(in) :: n1
complex(dp), pointer :: x(:)
type(c_ptr) :: p
p = fftw_alloc_complex(int(n1, c_size_t))
call c_f_pointer(p, x, [n1])
end function

function alloc3d(n1, n2, n3) result(x)
integer, intent(in) :: n1, n2, n3
complex(dp), pointer :: x(:, :, :)
type(c_ptr) :: p
p = fftw_alloc_complex(int(n1*n2*n3, c_size_t))
call c_f_pointer(p, x, [n1, n2, n3])
end function

subroutine free(x)
complex(dp), intent(in), target :: x(*)
call fftw_free(c_loc(x))
end subroutine

end module
