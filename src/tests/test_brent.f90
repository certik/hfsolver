program test_brent

use types, only: dp
use optimize, only: brent
use utils, only: assert
implicit none
real(dp) :: xmin, fxmin, eps
eps = 1e-12_dp

call brent(f1, -1._dp, 0.5_dp, 1._dp, eps, 20, xmin, fxmin)
call assert(abs(xmin - 0) < eps)
call brent(f1, -1._dp, 0._dp, 1._dp, eps, 20, xmin, fxmin)
call assert(abs(xmin - 0) < eps)

contains

real(dp) function f1(x) result(y)
real(dp), intent(in) :: x
y = x**2
end function

end program
