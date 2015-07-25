module mkl
use types, only: dp
implicit none
private
public fft3_inplace

interface

    integer(c_int) function fft3_inplace(x, n1, n2, n3) &
        bind(c, name="fft3_inplace")
    use iso_c_binding, only: c_int, c_double, c_double_complex
    implicit none
    integer(c_int), value, intent(in) :: n1, n2, n3
    complex(c_double_complex), intent(inout) :: x(n1*n2*n3)
    end function

end interface

end module
