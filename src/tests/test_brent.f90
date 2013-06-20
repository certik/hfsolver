program test_brent

use types, only: dp
use optimize, only: brent, bracket
use utils, only: assert
implicit none
real(dp) :: xmin, fxmin, eps, xa, xb, xc
eps = 1e-12_dp

call brent(f1, -1._dp, 0.5_dp, 1._dp, eps, 20, xmin, fxmin)
call assert(abs(xmin - 0) < eps)
call brent(f1, -1._dp, 0._dp, 1._dp, eps, 20, xmin, fxmin)
call assert(abs(xmin - 0) < eps)

call brent(f2, -4._dp, 0._dp, 4._dp, eps, 30, xmin, fxmin)
call assert(abs(xmin - 0.865474033101_dp) < eps)
call brent(f2, -4._dp, 1._dp, 4._dp, eps, 30, xmin, fxmin)
call assert(abs(xmin - 0.865474033101_dp) < eps)

xa = -4
xb = -3
call bracket(f2, xa, xb, xc, 100._dp, 20)
call brent(f2, xa, xb, xc, eps, 30, xmin, fxmin)
call assert(abs(xmin - 0.865474033101_dp) < eps)


contains

real(dp) function f1(x) result(y)
real(dp), intent(in) :: x
y = x**2
end function

real(dp) function f2(x) result(y)
real(dp), intent(in) :: x
y = (cos(x) - x**3)**2
end function

end program
