module petsc_
use iso_c_binding, only: c_int, c_double_complex
use types, only: dp
use utils, only: stop_error
implicit none
private
public petsc_init, petsc_finalize, solve

contains

integer(c_int) function solve(n, A, b, x) result(r)
integer(c_int), value, intent(in) :: n
complex(c_double_complex), intent(in) :: A(n*n), b(n)
complex(c_double_complex), intent(out) :: x(n)
r = 0
call stop_error("Not configured with PETSc.")
print *, n, A, b, x
end function

integer(c_int) function petsc_init() result(r)
r = 0
call stop_error("Not configured with PETSc.")
end function

integer(c_int) function petsc_finalize() result(r)
r = 0
call stop_error("Not configured with PETSc.")
end function

end module
