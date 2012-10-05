module c_hfsolver

use iso_c_binding, only: c_int, c_double
use special_functions, only: Inu_formula2, Knu_formula2
implicit none

contains

subroutine hfsolver_Inu(k, x, r) bind(c)
integer(c_int), intent(in) :: k
real(c_double), intent(in) :: x
real(c_double), intent(out) :: r
r = Inu_formula2(k, x)
end subroutine

subroutine hfsolver_Knu(k, x, r) bind(c)
integer(c_int), intent(in) :: k
real(c_double), intent(in) :: x
real(c_double), intent(out) :: r
r = Knu_formula2(k, x)
end subroutine

end module
