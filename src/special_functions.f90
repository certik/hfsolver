module special_functions
! Special functions
use types, only: dp
use constants, only: i_
use utils, only: stop_error
implicit none
private
public wigner3j, getgaunt, getgauntr

contains

real(dp) function fact(n) result(r)
integer, intent(in) :: n
r = gamma(n+1.0_dp)
end function

real(dp) function wigner3j(j1, j2, j3, m1, m2, m3) result(r)
! Calculates wigner3j symbols.
integer, intent(in) :: j1, j2, j3 ! angular momentum quantum numbers
integer, intent(in) :: m1, m2, m3 ! magnetic quantum numbers
integer :: k, k1, k2, l1, l2, l3, n1, n2, sgn
if ((abs(m1) > j1) .or. (abs(m2) > j2) .or. (abs(m3) > j3)) then
    print *, "Arguments:"
    print *, j1, j2, j3
    print *, m1, m2, m3
    call stop_error("wigner3j: invalid arguments")
end if
if ((j1 == 0) .and. (j2 == 0) .and. (j3 == 0)) then
    r = 1
    return
end if

l1 = j2 - j1 + j3
l2 = j1 - j2 + j3
l3 = j1 + j2 - j3
if ((m1 + m2 + m3 /= 0) .or. (l1 < 0) .or. (l2 < 0) .or. (l3 < 0)) then
    r = 0
    return
end if
n1 = j1 - m1
n2 = j2 + m2
k1 = max(0, j2-j3-m1, j1-j3+m2)
k2 = min(l3, n1, n2)
if (mod(k1 + j1 + j2 + m3, 2) /= 0) then
    sgn = -1
else
    sgn = 1
end if
r = 0
do k = k1, k2
    r = r + sgn * fr(l1, l1-n2+k) * fr(l2, l2-n1+k) * fr(l3, l3-k) / &
            (fact(k) * fact(n1-k) * fact(n2-k))
    sgn = -sgn
end do
r = r * sqrt(fr(j1+m1, l1) * fr(j2+m2, l2) * fr(j3+m3, l3) * &
        fr(j3-m3, 1+j1+j2+j3) * fact(j1-m1) * fact(j2-m2))
end function

real(dp) function fr(a, b) result(r)
! Returns a!/b!, for a,b >= 0
integer, intent(in) :: a, b
integer :: i
if (a < 0 .or. b < 0) call stop_error("fr(a,b): must satisfy a, b >= 0")
r = 1
if (a < b) then
    do i = a + 1, b
        r = r / i
    end do
else if (a > b) then
    do i = b + 1, a
        r = r * i
    end do
end if
end function

subroutine getgaunt(Lmax, ck)
! Calculates the Gaunt c^k(l1, m1, l2, m2) coefficients
! It returns an array of coefficients for all combinations of k, l1, m1, l2, m2:
!   ck(k, l1, m1, l2, m2) = c^k(l1, m1, l2, m2)
! Indices out of bounds mean that the coefficient is zero.
integer, intent(in) :: Lmax  ! max(l1, l2)
real(dp), allocatable, intent(out) :: ck(:, :, :, :, :) ! ck(k, l1, m1, l2, m2)
integer :: k, l1, m1, l2, m2
allocate(ck(0:2*Lmax, 0:Lmax, -Lmax:Lmax, 0:Lmax, -Lmax:Lmax))
ck = 0
do l1 = 0, Lmax
    do l2 = 0, Lmax
        do m1 = -l1, l1
            do m2 = -l2, l2
                do k = abs(l1-l2), l1+l2, 2
                    if (abs(m1-m2) > k) cycle
                    ck(k, l1, m1, l2, m2) = (-1)**(-m1) * &
                        sqrt(1._dp*(2*l1+1)*(2*l2+1)) * &
                        wigner3j(l1, k, l2, 0, 0, 0) * &
                        wigner3j(l1, k, l2, -m1, m1-m2, m2)
                end do
            end do
        end do
    end do
end do
end subroutine

subroutine getgauntr(Lmax, gr)
integer, intent(in) :: Lmax  ! max(l1, l2)
! gr(i, j, k, l, m, n) = <ij|kl|mn>_R
real(dp), allocatable, intent(out) :: gr(:, :, :, :, :, :)
real(dp), allocatable :: ck(:, :, :, :, :)
integer :: l1, m1, l2, m2, l3, m3
integer :: Lmax2
Lmax2 = 4*Lmax
allocate(ck(0:2*Lmax2, 0:Lmax2, -Lmax2:Lmax2, 0:Lmax2, -Lmax2:Lmax2))
call getgaunt(Lmax2, ck)
allocate(gr(0:Lmax, -Lmax:Lmax, 0:2*Lmax, -2*Lmax:2*Lmax, 0:Lmax, -Lmax:Lmax))
gr = 0
do l1 = 0, Lmax
    do l2 = 0, 2*Lmax
        do l3 = 0, Lmax
            do m1 = -l1, l1
                do m2 = -l2, l2
                    do m3 = -l3, l3
                        gr(l1, m1, l2, m2, l3, m3) = &
                            gauntr(l1, m1, l2, m2, l3, m3, ck, Lmax2)
                    end do
                end do
            end do
        end do
    end do
end do
end subroutine

real(dp) recursive function gauntr(l1, m1, l2, m2, l3, m3, ck, Lmax) result(r)
integer, intent(in) :: l1, m1, l2, m2, l3, m3, Lmax
real(dp), intent(in) :: ck(0:, 0:, -Lmax:, 0:, -Lmax:)
if (m1 /= 0 .and. m2 /= 0 .and. m3 /= 0) then
    ! Case A
    r = 2*ck(l2, l1, m2+m3, l3, m3) &
            * real(U(m1, m2+m3) * U(m2, m2) * U(m3, m3), dp) &
        + 2*ck(l2, l1, m2-m3, l3, -m3) &
            * real(U(m1, m2-m3) * U(m2, m2) * U(m3, -m3), dp)
else if (m1 /= 0 .and. m2 /= 0 .and. m3 == 0) then
    ! Case B
    r = 2*ck(l2, l1, m2, l3, 0) * real(U(m1, m2) * U(m2, m2), dp)
else if (m1 /= 0 .and. m2 == 0 .and. m3 /= 0) then
    ! Convert to case B
    r = gauntr(l1, m1, l3, m3, l2, m2, ck, Lmax)
else if (m1 == 0 .and. m2 /= 0 .and. m3 /= 0) then
    ! Convert to case B
    r = gauntr(l3, m3, l2, m2, l1, m1, ck, Lmax)
else if (m2 == 0 .and. m3 == 0) then
    ! Case C
    r = ck(l2, l1, 0, l3, 0) * delta(m1, 0)
else if (m1 == 0 .and. m3 == 0) then
    ! Convert to case C
    r = gauntr(l2, m2, l1, m1, l3, m3, ck, Lmax)
else if (m1 == 0 .and. m2 == 0) then
    ! Convert to case C
    r = gauntr(l3, m3, l1, m1, l2, m2, ck, Lmax)
else
    call stop_error("Internal error.")
end if

contains

integer function delta(a, b) result(r)
integer, intent(in) :: a, b
if (a == b) then
    r = 1
else
    r = 0
end if
end function

complex(dp) function U(mu, m) result(r)
integer, intent(in) :: mu, m
if (mu > 0) then
    r = (delta(mu, m) + (-1)**m * delta(mu, -m))/sqrt(2._dp)
else if (mu == 0) then
    r = delta(0, m)
else
    r = (delta(mu, -m) - (-1)**m * delta(mu, m))/(i_ * sqrt(2._dp))
end if
end function

end function

end module
