module integration

! Routines for calculating numerical integrals

use types, only: dp
use constants, only: pi
implicit none
private
public integrate_trapz_1, integrate_trapz_3, integrate_trapz_5, &
    integrate_trapz_7, integrate_simpson

contains

real(dp) pure function integrate_trapz_1(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)
real(dp) :: g(size(Rp))
integer :: N
N = size(Rp)
g = f * Rp
s = (g(1) + g(N)) / 2
s = s + sum(g(2:N-1))
end function

real(dp) pure function integrate_trapz_3(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)

real(dp) :: g(size(Rp))
integer :: N
N = size(Rp)
g = f * Rp
s = (9 * (g(1) + g(N)) + 28 * (g(2) + g(N-1)) + 23 * (g(3) + g(N-2))) / 24
s = s + sum(g(4:N-3))
end function

real(dp) pure function integrate_trapz_5(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)
real(dp) :: g(size(Rp))
integer :: N
N = size(Rp)
g = f * Rp
s = (  475 * (g(1) + g(N  )) &
    + 1902 * (g(2) + g(N-1)) &
    + 1104 * (g(3) + g(N-2)) &
    + 1586 * (g(4) + g(N-3)) &
    + 1413 * (g(5) + g(N-4)) &
    ) / 1440
s = s + sum(g(6:N-5))
end function

real(dp) pure function integrate_trapz_7(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)
real(dp) :: g(size(Rp))
integer :: N
N = size(Rp)
g = f * Rp
s = (  36799 * (g(1) + g(N  )) &
    + 176648 * (g(2) + g(N-1)) &
    +  54851 * (g(3) + g(N-2)) &
    + 177984 * (g(4) + g(N-3)) &
    +  89437 * (g(5) + g(N-4)) &
    + 130936 * (g(6) + g(N-5)) &
    + 119585 * (g(7) + g(N-6)) &
    ) / 120960
s = s + sum(g(8:N-7))
end function

real(dp) pure function integrate_simpson(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)
real(dp) :: g(size(Rp))
integer :: i, N
N = size(Rp)
g = f * Rp
s = 0
do i = 2, N-1, 2
    s = s + g(i-1) + 4*g(i) + g(i+1)
end do
s = s / 3
if (modulo(N, 2) == 0) then
    ! If N is even, add the last slice separately
    s = s + (5*g(N) + 8*g(N-1) - g(N-2)) / 12
end if
end function

end module
