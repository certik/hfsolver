program test_fm
use types, only: dp
use constants, only: pi
use special_functions, only: Fm
use utils, only: assert
implicit none

real(dp), allocatable :: F(:), G(:)
real(dp) :: x
real(dp), parameter :: eps = 1e-15_dp
real(dp), parameter :: xlist_high(*) = [1e-6_dp, 0.1_dp, 0.5_dp, 2._dp, &
    10._dp, 15._dp, 19._dp, 200._dp, 500._dp]
real(dp), parameter :: xlist_low(*) = [20._dp, 21._dp, 25._dp, 30._dp, &
    50._dp, 100._dp, 150._dp]
integer :: m, maxm, i
! We use downward recursion with initial value 0, so maxm must be very high to
! get accurate results:
maxm = 1000
allocate(F(0:maxm), G(0:maxm))

do i = 1, size(xlist_high)
    x = xlist_high(i)
    call downwards(maxm, x, F)
    call assert(abs(F(0) - 0.5_dp * sqrt(pi/x) * erf(sqrt(x))) < eps)
    ! With these "x" we can go all the way to m=500 with high accuracy:
    do m = 0, 500
        call Fm(m, x, G(:m))
        call assert(all(abs(G(:m) - F(:m)) < eps))
    end do
end do

do i = 1, size(xlist_low)
    x = xlist_low(i)
    call downwards(maxm, x, F)
    call assert(abs(F(0) - 0.5_dp * sqrt(pi/x) * erf(sqrt(x))) < eps)
    ! With these "x" we can only go to m=50 with high accuracy:
    do m = 0, 50
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

end program
