module mkl
use types, only: dp
use utils, only: stop_error
implicit none
private
public fft3_inplace, ifft3_inplace

interface

    integer(c_int) function mkl_fft3_inplace(x, n1, n2, n3) &
        bind(c, name="mkl_fft3_inplace")
    use iso_c_binding, only: c_int, c_double, c_double_complex
    implicit none
    integer(c_int), value, intent(in) :: n1, n2, n3
    complex(c_double_complex), intent(inout) :: x(n1*n2*n3)
    end function

    integer(c_int) function mkl_ifft3_inplace(x, n1, n2, n3) &
        bind(c, name="mkl_ifft3_inplace")
    use iso_c_binding, only: c_int, c_double, c_double_complex
    implicit none
    integer(c_int), value, intent(in) :: n1, n2, n3
    complex(c_double_complex), intent(inout) :: x(n1*n2*n3)
    end function

end interface

contains

subroutine fft3_inplace(x)
complex(dp), intent(inout) :: x(:, :, :)
integer :: r
r = mkl_fft3_inplace(x, size(x, 1), size(x, 2), size(x, 3))
if (r /= 0) call stop_error("mkl_fft3_inplace returned an error.")
end subroutine

subroutine ifft3_inplace(x)
complex(dp), intent(inout) :: x(:, :, :)
integer :: r
r = mkl_ifft3_inplace(x, size(x, 1), size(x, 2), size(x, 3))
if (r /= 0) call stop_error("mkl_ifft3_inplace returned an error.")
end subroutine

end module
