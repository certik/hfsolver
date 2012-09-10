program test_gaunt

! Compare the output from this test with the gaunt_log file, which was manually
! checked against the table 1^6, page 178, 179 in Condon & Shortley, The Theory
! of Atomic Spectra --- all numbers there are correct.

use types, only: dp
use special_functions, only: getgaunt
use utils, only: assert
use hfprint, only: lmap
implicit none

integer :: Lmax
real(dp), allocatable :: ck(:, :, :, :, :)
integer :: l, lp, m, mp, k

Lmax = 3
call getgaunt(Lmax, ck)
print *, "l + l' odd"
do l = 0, Lmax
    do lp = l, Lmax
        if (mod(l + lp, 2) == 0) cycle
        print "(a, 3i9)", "Dk:       ", Dk(1, l, lp), Dk(3, l, lp), Dk(5, l, lp)
        do m = l, 0, -1
            do mp = lp, 0, -1
                print "(2a2,2i3,3i9)", lmap(l), lmap(lp), m, mp, &
                    f(1, ck(1, l, m, lp, mp)), &
                    f(3, ck(3, l, m, lp, mp)), &
                    f(5, ck(5, l, m, lp, mp))
                do k = 0, 2*Lmax
                    call assert(ck(k, l, m, lp, mp)-ck(k, l, -m, lp, -mp) &
                        < 1e-10_dp)
                end do
            end do
        end do
        do m = l, 1, -1
            do mp = lp, 1, -1
                print "(2a2,2i3,3i9)", lmap(l), lmap(lp), m, -mp, &
                    f(1, ck(1, l, m, lp, -mp)), &
                    f(3, ck(3, l, m, lp, -mp)), &
                    f(5, ck(5, l, m, lp, -mp))
                do k = 0, 2*Lmax
                    call assert(ck(k, l, m, lp, -mp)-ck(k, l, -m, lp, mp) &
                        < 1e-10_dp)
                end do
            end do
        end do
    end do
end do
print *
print *, "l + l' even"
do l = 0, Lmax
    do lp = l, Lmax
        if (mod(l + lp, 2) == 1) cycle
        print "(a, 4i10)", "Dk:       ", Dk(0, l, lp), Dk(2, l, lp), &
            Dk(4, l, lp), Dk(6, l, lp)
        do m = l, 0, -1
            do mp = lp, 0, -1
                print "(2a2,2i3,4i10)", lmap(l), lmap(lp), m, mp, &
                    f(0, ck(0, l, m, lp, mp)), &
                    f(2, ck(2, l, m, lp, mp)), &
                    f(4, ck(4, l, m, lp, mp)), &
                    f(6, ck(6, l, m, lp, mp))
                do k = 0, 2*Lmax
                    call assert(ck(k, l, m, lp, mp)-ck(k, l, -m, lp, -mp) &
                        < 1e-10_dp)
                end do
            end do
        end do
        do m = l, 1, -1
            do mp = lp, 1, -1
                print "(2a2,2i3,4i10)", lmap(l), lmap(lp), m, -mp, &
                    f(0, ck(0, l, m, lp, -mp)), &
                    f(2, ck(2, l, m, lp, -mp)), &
                    f(4, ck(4, l, m, lp, -mp)), &
                    f(6, ck(6, l, m, lp, -mp))
                do k = 0, 2*Lmax
                    call assert(ck(k, l, m, lp, -mp)-ck(k, l, -m, lp, mp) &
                        < 1e-10_dp)
                end do
            end do
        end do
    end do
end do

contains

integer function Dk(k, l, lp)
integer, intent(in) :: k, l, lp
Dk = 0
if (l == 0 .and. lp == 0) then
    if (k == 0) Dk = 1
else if (l == 0 .and. lp == 1) then
    if (k == 1) Dk = 3
else if (l == 0 .and. lp == 2) then
    if (k == 2) Dk = 5
else if (l == 0 .and. lp == 3) then
    if (k == 3) Dk = 7
else if (l == 1 .and. lp == 1) then
    if (k == 0) then
        Dk = 1
    else if (k == 2) then
        Dk = 25
    end if
else if (l == 1 .and. lp == 2) then
    if (k == 1) then
        Dk = 15
    else if (k == 3) then
        Dk = 245
    end if
else if (l == 1 .and. lp == 3) then
    if (k == 2) then
        Dk = 175
    else if (k == 4) then
        Dk = 189
    end if
else if (l == 2 .and. lp == 2) then
    if (k == 0) then
        Dk = 1
    else if (k == 2) then
        Dk = 49
    else if (k == 4) then
        Dk = 441
    end if
else if (l == 2 .and. lp == 3) then
    if (k == 1) then
        Dk = 35
    else if (k == 3) then
        Dk = 315
    else if (k == 5) then
        Dk = 15246
    end if
else if (l == 3 .and. lp == 3) then
    if (k == 0) then
        Dk = 1
    else if (k == 2) then
        Dk = 225
    else if (k == 4) then
        Dk = 1089
    else if (k == 6) then
        Dk = 736164
    end if
end if
end function

integer function f(k, a) result(r)
integer, intent(in) :: k
real(dp), intent(in) :: a
real(dp) :: rr
if (a == 0) then
    r = 0
else
    rr = abs(a) / a  * a**2 * Dk(k, l, lp)
    call assert(Dk(k, l, lp) /= 0)
    r = round(rr)
    ! This tests that the calculated value is 1e-10 from the exact integer:
    call assert(abs(r - rr) < 1e-10_dp)
end if
end function

integer pure function round(x) result(i)
real(dp), intent(in) :: x
if (x == 0) then
    i = 0
else
    i = int(x + abs(x)/x * 0.5)
end if
end function

end program
