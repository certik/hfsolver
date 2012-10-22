module petsc_
use types, only: dp
implicit none
private
public petsc_init, petsc_finalize, solve

interface

    integer(c_int) function solve(n, A, b, x) bind(c, name="petsc_solve")
    use iso_c_binding, only: c_int, c_double, c_double_complex
    implicit none
    integer(c_int), value, intent(in) :: n
    complex(c_double_complex), intent(in) :: A(n*n), b(n)
    complex(c_double_complex), intent(out) :: x(n)
    end function

    integer(c_int) function petsc_init() bind(c)
    use iso_c_binding, only: c_int
    implicit none
    end function

    integer(c_int) function petsc_finalize() bind(c)
    use iso_c_binding, only: c_int
    implicit none
    end function

end interface

end module
