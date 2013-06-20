program test_brent

use types, only: dp
use optimize, only: brent
implicit none
real(dp) :: xmin, fxmin

call brent(f1, -1._dp, 0._dp, 1._dp, 1e-12_dp, 200, xmin, fxmin)
print *, xmin, fxmin

contains

real(dp) function f1(x) result(y)
real(dp), intent(in) :: x
y = x**2
end function

end program
