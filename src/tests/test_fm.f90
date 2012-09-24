program test_fm
use types, only: dp
use constants, only: pi
use special_functions, only: Fm
use utils, only: assert
implicit none

real(dp), allocatable :: F(:), G(:)
real(dp) :: x
real(dp), parameter :: eps = 1e-15_dp
real(dp), parameter :: eps_rel = 1e-14_dp
real(dp), parameter :: xlist(*) = [1e-6_dp, 0.1_dp, 0.5_dp, 2._dp, 10._dp, &
    15._dp, 19._dp, 20._dp, 21._dp, 25._dp, 30._dp, 50._dp, 100._dp, 150._dp, &
    200._dp, 500._dp, 2000._dp, 20000._dp, 2e5_dp]
integer :: m, maxm, i
! We use downward recursion with initial value 0, so maxm must be very high to
! get accurate results:
maxm = 1000
allocate(F(0:maxm), G(0:maxm))

do i = 1, size(xlist)
    x = xlist(i)
    if (x < 510) then
        call downwards(maxm, x, F)
    else
        call asympt(maxm, x, F)
    end if
    call assert(abs(F(0) - 0.5_dp * sqrt(pi/x) * erf(sqrt(x))) < eps)
    call assert(abs(F(0) - 0.5_dp * sqrt(pi/x) * erf(sqrt(x))) / F(0) < eps_rel)
    do m = 0, 500
        call Fm(m, x, G(:m))
        call assert(all(abs(G(:m) - F(:m)) < eps))
    end do
end do

contains

subroutine downwards(maxm, t, F)
integer, intent(in) :: maxm
real(dp), intent(in) :: t
real(dp), intent(out) :: F(0:)
integer :: m
F(maxm) = 0
do m = maxm-1, 0, -1
    F(m) = (2*t*F(m + 1) + exp(-t)) / (2*m + 1)
end do
end subroutine

subroutine asympt(maxm, t, F)
integer, intent(in) :: maxm
real(dp), intent(in) :: t
real(dp), intent(out) :: F(0:)
integer :: m
do m = 0, maxm
    ! This can overflow (NaN):
    !F(m) = gamma(m+0.5_dp) / (2*t**(m+0.5_dp))
    ! so this is equivalent, but works for large "t" and "m":
    F(m) = exp(log_gamma(m+0.5_dp) - (m+0.5_dp)*log(t)) / 2
end do
end subroutine

end program
