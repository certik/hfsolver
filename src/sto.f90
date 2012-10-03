module sto
use types, only: dp
use constants, only: pi
use utils, only: stop_error
use scf, only: ijkl2intindex
use utils, only: str, assert
use hartree_screening, only: hartree_y
use special_functions, only: wigner3j, getgaunt
use utils, only: loadtxt
use openmp, only: omp_get_thread_num
use quadrature, only: gauss_pts, gauss_wts
!use debye, only: Vk
!use debye, only: Vk => Sk
!use debye, only: Vk => Vk2
use debye, only: Vk => Vk3
implicit none
private
public get_basis, get_basis2, stoints, stoints2, get_values, slater_sto_gauss, &
    slater_sto_screen, sto_V_screen

! Array of factorials: fact(n) = n!
! You have to call calc_factorials() to initialize it first.
real(dp), allocatable :: fact(:)

! only used in spec()
real(dp), allocatable :: xiq(:), wtq(:)
! only used in gauss_ab()
integer, parameter :: Nq = 64
real(dp), allocatable :: xiq3(:), wtq3(:)

contains

subroutine get_basis(nl, zl, ll, nlist, zetalist, llist, mlist)
integer, intent(in) :: nl(:), ll(:)
real(dp), intent(in) :: zl(:)
integer, allocatable, intent(out) :: nlist(:), llist(:), mlist(:)
real(dp), allocatable, intent(out) :: zetalist(:)
integer :: i, j, n, l, m, a
real(dp) :: zeta
a = sum(2*ll+1)
allocate(nlist(a), zetalist(a), llist(a), mlist(a))
i = 1
j = 1
do j = 1, size(nl)
    n = nl(j)
    zeta = zl(j)
    l = ll(j)
    do m = -l, l
        nlist(i) = n
        zetalist(i) = zeta
        llist(i) = l
        mlist(i) = m
        i = i + 1
    end do
end do
end subroutine

subroutine calc_factorials(maxn, fact)
! Returns fact(n)=n!  for n = 0..maxn
integer, intent(in) :: maxn
real(dp), intent(out), allocatable :: fact(:)
integer :: n
allocate(fact(0:maxn))
do n = 0, maxn
    fact(n) = gamma(n + 1._dp)
end do
end subroutine

real(dp) pure function sto_norm(n, zeta) result(r)
integer, intent(in) :: n
real(dp), intent(in) :: zeta
r = sqrt((2*zeta)**(2*n+1) / fact(2*n))
end function

real(dp) function overlap(n, zeta, i, j) result(r)
integer, intent(in) :: n(:)
real(dp), intent(in) :: zeta(:)
integer, intent(in) :: i, j
r = sto_norm(n(i), zeta(i)) * sto_norm(n(j), zeta(j)) * &
        fact(n(i)+n(j)) / (zeta(i)+zeta(j))**(n(i)+n(j)+1)
end function

real(dp) function potential(n, zeta, Z, i, j) result(r)
integer, intent(in) :: n(:), Z
real(dp), intent(in) :: zeta(:)
integer, intent(in) :: i, j
r = -Z * sto_norm(n(i), zeta(i)) * sto_norm(n(j), zeta(j)) * &
        fact(n(i)+n(j)-1) / (zeta(i)+zeta(j))**(n(i)+n(j))
end function

real(dp) function potential_screen(n, zeta, Z, i, j, D) result(r)
integer, intent(in) :: n(:), Z
real(dp), intent(in) :: zeta(:)
integer, intent(in) :: i, j
real(dp), intent(in) :: D
r = -Z * sto_norm(n(i), zeta(i)) * sto_norm(n(j), zeta(j)) * &
        fact(n(i)+n(j)-1) / (zeta(i)+zeta(j)+1/D)**(n(i)+n(j))
end function

real(dp) function kinetic(n, zeta, l, i, j) result(r)
integer, intent(in) :: n(:)
real(dp), intent(in) :: zeta(:)
integer, intent(in) :: l
integer, intent(in) :: i, j
r = 0.5_dp * sto_norm(n(i), zeta(i)) * sto_norm(n(j), zeta(j)) * ( &
        (n(i)*n(j) + l*(l+1)) * fact(n(i)+n(j)-2) / &
            (zeta(i)+zeta(j))**(n(i)+n(j)-1) &
        - (n(i)*zeta(j)+n(j)*zeta(i)) * fact(n(i)+n(j)-1) / &
            (zeta(i)+zeta(j))**(n(i)+n(j)) &
        + zeta(i)*zeta(j) * fact(n(i)+n(j)) / &
            (zeta(i)+zeta(j))**(n(i)+n(j)+1) )
end function

real(dp) pure function H(n, zeta, kk, i, j, k, l) result(r)
integer, intent(in) :: n(:)
real(dp), intent(in) :: zeta(:)
integer, intent(in) :: kk, i, j, k, l
integer :: nu
r = 0
do nu = 0, n(i)+n(k)-kk-1
    r = r + fact(n(j)+n(l)+kk+nu) * (zeta(i)+zeta(k))**nu &
        / (fact(nu) * (zeta(i)+zeta(j)+zeta(k)+zeta(l))**(n(j)+n(l)+kk+nu+1))
end do
end function

real(dp) pure function slater(n, zeta, kk, i, j, k, l) result(r)
! Calculates slater integral R^kk(i, j, k, l) of the STO orbitals given in the
! 'n' and 'zeta' arrays. The indices <ij|kl> are in the physical notation.
integer, intent(in) :: n(:)
real(dp), intent(in) :: zeta(:)
integer, intent(in) :: kk, i, j, k, l
if (n(i)+n(k)-kk-1 < 0 .or. n(j)+n(l)-kk-1 < 0) then
    ! These integrals are never needed due to the winger 3j selection rule, so
    ! we don't calculate them and simply return zero. (Note that the formulas
    ! that we use wouldn't work anyway for this case, so one would have to use
    ! a different method, but that's a different issue.)
    r = 0
else
    r = sto_norm(n(i), zeta(i)) * sto_norm(n(j), zeta(j)) * &
            sto_norm(n(k), zeta(k)) * sto_norm(n(l), zeta(l)) * ( &
                fact(n(i)+n(k)-kk-1) / (zeta(i)+zeta(k))**(n(i)+n(k)-kk) &
                    * H(n, zeta, kk, i, j, k, l) &
              + fact(n(j)+n(l)-kk-1) / (zeta(j)+zeta(l))**(n(j)+n(l)-kk) &
                    * H(n, zeta, kk, j, i, l, k))
end if
end function

subroutine stoints(Z, nlist, zetalist, llist, mlist, S, T, V, int2)
integer, intent(in) :: Z
integer, intent(in) :: nlist(:), llist(:), mlist(:)
real(dp), intent(in) :: zetalist(:)
real(dp), intent(out) :: S(:, :), T(:, :), V(:, :)
real(dp), intent(out) :: int2(:)
integer :: i, j, k, l, Lmax, k_
integer :: l1, l1p, l2, l2p, m1, m1p, m2, m2p
real(dp) :: r, fac
real(dp), allocatable :: ck(:, :, :, :, :)
! Precalculate factorials. The maximum factorial needed is given by 4*maxn-1:
call calc_factorials(maxval(nlist)*4-1, fact)
! Precalculate the radial integrals
Lmax = maxval(llist)
! Construct S, T, V
S = 0
T = 0
V = 0
print *, "Calculating S, T, V, ...."
do i = 1, size(nlist)
    do j = 1, size(nlist)
        if (llist(i) == llist(j) .and. mlist(i) == mlist(j)) then
            S(i, j) = overlap(nlist, zetalist, i, j)
            T(i, j) = kinetic(nlist, zetalist, llist(i), i, j)
            V(i, j) = potential(nlist, zetalist, Z, i, j)
        end if
    end do
end do
print *, "Calculating Gaunt coefficients..."
call getgaunt(Lmax, ck)
!call print_gaunt(Lmax, ck)
print *, "Calculating int2..."
int2 = 0
do i = 1, size(nlist)
    print *, 100._dp * i / size(nlist), "%"
    do j = i, size(nlist)
        do k = 1, size(nlist)
            do l = k, size(nlist)
                if ((i-1)*i/2 + j < (k-1)*k/2 + l) cycle
                l1 = llist(i); l1p = llist(j)
                l2 = llist(k); l2p = llist(l)
                m1 = mlist(i); m1p = mlist(j)
                m2 = mlist(k); m2p = mlist(l)
                if (m1+m2 /= m1p+m2p) cycle
!                m1 = 0
!                m2 = 0
!                m1p = 0
!                m2p = 0
                r = 0
                do k_ = max(abs(l1-l1p), abs(l2-l2p)), min(l1+l1p, l2+l2p)
                    fac = ck(k_, l1, m1, l1p, m1p) * ck(k_, l2p, m2p, l2, m2)
                    if (fac == 0) cycle
                    r = r + fac * slater(nlist, zetalist, k_, i, k, j, l)
                end do
                int2(ijkl2intindex(i, j, k, l)) = r
            end do
        end do
    end do
end do
end subroutine

subroutine stoints2(Z, nbfl, nlist, zetalist, S, T, V, slater_)
! Just like stoints, but it only handles the radial part.
! To get the full 3D integrals, use slater2int22 and radialC2C.
integer, intent(in) :: Z
integer, intent(in) :: nbfl(0:), nlist(:, 0:)
real(dp), intent(in) :: zetalist(:, 0:)
real(dp), allocatable, intent(out) :: S(:, :, :), T(:, :, :), V(:, :, :), &
    slater_(:, :)
integer,  dimension(sum(nbfl)) :: nl
real(dp), dimension(sum(nbfl)) :: zl
integer :: m, n, mu, nu, i, j, k, l, Lmax, k_, ijkl, ndof
integer, allocatable :: intindex(:, :, :, :)
! Precalculate factorials. The maximum factorial needed is given by 4*maxn-1:
call calc_factorials(maxval(nlist)*4-1, fact)
! Precalculate the radial integrals
Lmax = ubound(nbfl, 1)
ndof = sum(nbfl)
n = maxval(nbfl)
m = ndof*(ndof+1)/2
allocate(S(n, n, 0:Lmax), T(n, n, 0:Lmax), V(n, n, 0:Lmax))
allocate(slater_(m*(m+1)/2, 0:2*Lmax))
! Construct S, T, V
S = 0
T = 0
V = 0
print *, "Calculating S, T, V, ...."
do l = 0, ubound(nbfl, 1)
    do mu = 1, nbfl(l)
        do nu = 1, nbfl(l)
            S(mu, nu, l) = overlap  (nlist(:, l), zetalist(:, l),    mu, nu)
            T(mu, nu, l) = kinetic  (nlist(:, l), zetalist(:, l), l, mu, nu)
            V(mu, nu, l) = potential(nlist(:, l), zetalist(:, l), Z, mu, nu)
        end do
    end do
end do
print *, "Calculating Slater integrals..."
slater_ = 0
j = 1
do l = 0, ubound(nbfl, 1)
    nl(j:j+nbfl(l)-1) = nlist   (:nbfl(l), l)
    zl(j:j+nbfl(l)-1) = zetalist(:nbfl(l), l)
    j = j + nbfl(l)
end do
n = size(nl)
allocate(intindex(n, n, n, n))
call create_index(intindex)
!$omp parallel default(none) shared(n, nbfl, nl, zl, slater_, intindex) private(i, j, k, l, k_, ijkl)
!$omp do schedule(dynamic)
do i = 1, n
    do j = i, n
        do k = 1, n
            do l = k, min((i-1)*i/2 + j - (k-1)*k/2, n)
                ijkl = intindex(l, k, j, i)
                do k_ = 0, 2*ubound(nbfl, 1)
                    slater_(ijkl, k_) = slater(nl, zl, k_, i, k, j, l)
                end do
            end do
        end do
    end do
    ! char(13) is carriage return, so we keep overwriting the percentage
!    write (*, "(a1, f5.1,'%')", advance="no") char(13), &
!        100._dp * ijkl / size(slater_, 1)
end do
!$omp end do
!$omp end parallel
!print *
end subroutine

subroutine sto_V_screen(Z, nbfl, nlist, zetalist, V, D)
! Calculates the potential matrix V with Debye screening
integer, intent(in) :: Z
integer, intent(in) :: nbfl(0:), nlist(:, 0:)
real(dp), intent(in) :: zetalist(:, 0:)
real(dp), allocatable, intent(out) :: V(:, :, :)
real(dp), intent(in) :: D ! Debye screening parameter
integer ::  n, l, mu, nu, Lmax
! Precalculate factorials. The maximum factorial needed is given by 4*maxn-1:
call calc_factorials(maxval(nlist)*4-1, fact)
! Precalculate the radial integrals
Lmax = ubound(nbfl, 1)
n = maxval(nbfl)
allocate(V(n, n, 0:Lmax))
V = 0
print *, "Calculating V with screening ...."
do l = 0, ubound(nbfl, 1)
    do mu = 1, nbfl(l)
        do nu = 1, nbfl(l)
            V(mu, nu, l) = potential_screen(nlist(:, l), zetalist(:, l), &
                Z, mu, nu, D)
        end do
    end do
end do
end subroutine

real(dp) function factor(l1, l1p, l2, l2p, m1, m1p, m2, m2p, k) result(r)
integer, intent(in) :: l1, l1p, l2, l2p, m1, m1p, m2, m2p, k
integer :: q
r = 0
do q = -k, k
    r = r + (-1)**(m1+m2+q) &
        * wigner3j(l1, l1p, k, -m1, m1p, +q) &
        * wigner3j(l2, l2p, k, -m2, m2p, -q)
end do
if (abs(r-factor2(l1, l1p, l2, l2p, m1, m1p, m2, m2p, k)) > 1e-9_dp) then
    stop "FAIL"
end if
end function

real(dp) function factor2(l1, l1p, l2, l2p, m1, m1p, m2, m2p, k) result(r)
integer, intent(in) :: l1, l1p, l2, l2p, m1, m1p, m2, m2p, k
if (k >= abs(m1-m1p) .and. m1+m2-m1p-m2p == 0) then
    r = (-1)**(m1+m2p) &
        * wigner3j(l1, l1p, k, -m1, m1p, m1-m1p) &
        * wigner3j(l2, l2p, k, -m2, m2p, m2-m2p)
else
    r = 0
end if
end function

real(dp) function factor3(l1, l1p, l2, l2p, m1, m1p, m2, m2p, k, Lmax, &
        wigner3j__) result(r)
integer, intent(in) :: l1, l1p, l2, l2p, m1, m1p, m2, m2p, k, Lmax
real(dp), intent(in) :: wigner3j__(0:, 0:, 0:, -Lmax:, -Lmax:, -2*Lmax:)
if (k >= abs(m1-m1p) .and. m1+m2-m1p-m2p == 0) then
    r = (-1)**(m1+m2p) &
        * wigner3j__(l1, l1p, k, -m1, m1p, m1-m1p) &
        * wigner3j__(l2, l2p, k, -m2, m2p, m2-m2p)
else
    r = 0
end if
end function

subroutine print_gaunt(Lmax, ck)
integer, intent(in) :: Lmax
real(dp), intent(in) :: ck(0:, 0:, -Lmax:, 0:, -Lmax:)
integer :: i, j, k, m1, m2
do k = lbound(ck, 1), ubound(ck, 1)
    do i = lbound(ck, 2), ubound(ck, 2)
        do m1 = lbound(ck, 3), ubound(ck, 3)
            do j = lbound(ck, 4), ubound(ck, 4)
                do m2 = lbound(ck, 5), ubound(ck, 5)
                    if (ck(k, i, m1, j, m2) == 0) cycle
                    print "('   ', i3, i3, i3, i3, i3, es24.16)", k, i, &
                        m1, j, m2, ck(k, i, m1, j, m2)
                end do
            end do
        end do
    end do
end do
end subroutine

function get_values(nbfl, nl, zl, C, R) result(u)
integer, intent(in) :: nl(:, 0:), nbfl(0:)
real(dp), intent(in) :: zl(:, 0:), C(:, :, 0:), R(:)
! u(:, n, l) are the values of eigenvector (n, l)
real(dp) :: u(size(R), size(C, 2), 0:ubound(C, 3))
integer :: n, l
call calc_factorials(maxval(nl)*2, fact)
u = 0
do l = 0, ubound(nbfl, 1)
    do n = 1, nbfl(l)
        u(:, n, l) = sol_eval(nl(:nbfl(l), l), zl(:nbfl(l), l), &
                C(:nbfl(l), n, l), R)
    end do
end do
end function

function sol_eval(n, z, C, r) result(u)
integer, intent(in) :: n(:)
real(dp), intent(in) :: z(:), C(:), r(:)
real(dp) :: u(size(R))
integer :: i
u = 0
do i = 1, size(n)
    u = u + C(i) * sqrt((2*z(i))**(2*n(i)+1) / fact(2*n(i))) * r**n(i) * &
            exp(-z(i)*r)
end do
end function


subroutine get_basis2(Z, nbfl, nl, zl, focc)
integer, intent(in) :: Z
integer, intent(out), allocatable :: nl(:, :), nbfl(:)
real(dp), intent(out), allocatable :: zl(:, :), focc(:, :)

if (Z == 2) then
    allocate(nl(5, 0:0), zl(5, 0:0), nbfl(0:0), focc(5, 0:0))
    focc = 0
    nbfl = 0
    focc(1, 0) = 2
    nbfl(0) = 5
    nl(:5, 0) = [1, 1, 1, 1, 1]
    zl(:5, 0) = [1.41714_dp, 2.37682_dp, 4.39628_dp, 6.52699_dp, 7.94252_dp]
else if (Z == 4) then
    allocate(nl(6, 0:0), zl(6, 0:0), nbfl(0:0), focc(6, 0:0))
    focc = 0
    nbfl = 0
    focc(:2, 0) = [2, 2]
    nbfl(0) = 6
    nl(:6, 0) = [1, 1, 2, 2, 2, 2]
    zl(:6, 0) = [3.47116_dp, 6.36861_dp, 0.77820_dp, 0.94067_dp, 1.48725_dp, &
        2.71830_dp]
else if (Z == 10) then
    allocate(nl(6, 0:1), zl(6, 0:1), nbfl(0:1), focc(6, 0:1))
    focc = 0
    nl = 0
    zl = 0
    ! s
    focc(:2, 0) = [2, 2]
    nbfl(0) = 6
    nl(:6, 0) = [1, 1, 2, 2, 2, 2]
    zl(:6, 0) = [9.48486_dp, 15.56590_dp, 1.96184_dp, 2.86423_dp, 4.82530_dp, &
        7.79242_dp]

    ! p
    focc(:1, 1) = [6]
    nbfl(1) = 4
    nl(:4, 1) = [2, 2, 2, 2]
    zl(:4, 1) = [1.45208_dp, 2.38168_dp, 4.48489_dp, 9.13464_dp]
else if (Z == 12) then
    allocate(nl(8, 0:1), zl(8, 0:1), nbfl(0:1), focc(8, 0:1))
    focc = 0
    nl = 0
    zl = 0
    ! s
    focc(:3, 0) = [2, 2, 2]
    nbfl(0) = 8
    nl(:8, 0) = [1, 3, 3, 3, 3, 3, 3, 3]
    zl(:8, 0) = [12.01140_dp, 13.91620_dp, 9.48612_dp, 6.72188_dp, 4.24466_dp, &
        2.53466_dp, 1.46920_dp, 0.89084_dp]

    ! p
    focc(:1, 1) = [6]
    nbfl(1) = 5
    nl(:5, 1) = [2, 4, 4, 4, 4]
    zl(:5, 1) = [5.92580_dp, 7.98979_dp, 5.32964_dp, 3.71678_dp, 2.59986_dp]
else if (Z == 18) then
    allocate(nl(8, 0:1), zl(8, 0:1), nbfl(0:1), focc(8, 0:1))
    focc = 0
    nl = 0
    zl = 0
    ! s
    focc(:3, 0) = [2, 2, 2]
    nbfl(0) = 8
    nl(:8, 0) = [1, 3, 3, 3, 3, 3, 3, 3]
    zl(:8, 0) = [18.01640_dp, 22.04650_dp, 16.08250_dp, 11.63570_dp, &
        7.70365_dp, 4.87338_dp, 3.32987_dp, 2.02791_dp]

    ! p
    focc(:2, 1) = [6, 6]
    nbfl(1) = 8
    nl(:8, 1) = [2, 4, 4, 4, 4, 4, 4, 4]
    zl(:8, 1) = [9.05477_dp, 15.54410_dp, 12.39970_dp, 8.56120_dp, &
        5.94658_dp, 3.42459_dp, 1.96709_dp, 1.06717_dp]
else if (Z == 54) then
    allocate(nl(11, 0:2), zl(11, 0:2), nbfl(0:2), focc(5, 0:2))
    focc = 0
    nl = 0
    zl = 0
    ! s
    focc(:5, 0) = [2, 2, 2, 2, 2]
    nbfl(0) = 11
    nl(:11, 0) = [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 5]
    zl(:11, 0) = [55.30720_dp, 37.80730_dp, 27.92970_dp, 23.69210_dp, &
        15.03530_dp, 12.67230_dp, 7.60195_dp, 5.73799_dp, 4.17583_dp, &
        2.99772_dp, 1.98532_dp]

    ! p
    focc(:4, 1) = [6, 6, 6, 6]
    nbfl(1) = 9
    nl(:9, 1) = [2, 2, 3, 3, 4, 4, 5, 5, 5]
    zl(:9, 1) = [34.88440_dp, 23.30470_dp, 12.54120_dp, 12.02300_dp, &
        7.72390_dp, 5.40562_dp, 3.32661_dp, 2.09341_dp, 1.36686_dp]

    ! d
    focc(:2, 2) = [10, 10]
    nbfl(2) = 5
    nl(:5, 2) = [3, 3, 4, 4, 4]
    zl(:5, 2) = [20.08240_dp, 11.78600_dp, 7.30842_dp, 4.88400_dp, 3.19850_dp]
else
    call stop_error("Z=" // str(Z) // " not implemented")
end if
end subroutine


subroutine slater_sto_gauss(nbfl, nlist, zetalist, slater_)
! Just like stoints2, but only the slater integral and uses the Gaussian
! integration to do it.
integer, intent(in) :: nbfl(0:), nlist(:, 0:)
real(dp), intent(in) :: zetalist(:, 0:)
real(dp), allocatable, intent(out) :: slater_(:, :)
integer,  dimension(sum(nbfl)) :: nl
real(dp), dimension(sum(nbfl)) :: zl
integer :: m, n, i, j, k, l, k_, ijkl, ndof, Lmax
real(dp), allocatable :: d(:, :)
! Precalculate factorials. The maximum factorial needed is given by 4*maxn-1:
call calc_factorials(maxval(nlist)*4-1, fact)
if (.not. allocated(xiq)) then
    !call loadtxt("lag06.txt", d)
    !call loadtxt("lag20.txt", d)
    call loadtxt("lag52.txt", d)
    allocate(xiq(size(d, 1)), wtq(size(d, 1)))
    xiq = d(:, 1)
    wtq = d(:, 2)
end if
if (.not. allocated(xiq3)) then
    allocate(xiq3(Nq), wtq3(Nq))
    xiq3 = gauss_pts(Nq)
    wtq3 = gauss_wts(Nq)
end if
print *, "Calculating Slater integrals..."
Lmax = ubound(nbfl, 1)
ndof = sum(nbfl)
m = ndof*(ndof+1)/2
allocate(slater_(m*(m+1)/2, 0:2*Lmax))
slater_ = 0
j = 1
do l = 0, ubound(nbfl, 1)
    nl(j:j+nbfl(l)-1) = nlist   (:nbfl(l), l)
    zl(j:j+nbfl(l)-1) = zetalist(:nbfl(l), l)
    j = j + nbfl(l)
end do
n = size(nl)
ijkl = 0
do i = 1, n
    do j = 1, i
        do k = 1, n
            do l = 1, k
                if ((i-1)*i/2 + j < (k-1)*k/2 + l) cycle
                ijkl = ijkl + 1
                ! This check can be used to verify, that the "ijkl" index is
                ! calculated properly:
                !call assert(ijkl2intindex(i, j, k, l) == ijkl)
                do k_ = 0, 2*ubound(nbfl, 1)
                    !slater_(ijkl, k_) = slater(nl, zl, k_, i, k, j, l)
                    slater_(ijkl, k_) = slater_gauss(nl, zl, k_, i, k, j, l) &
                        + slater_gauss(nl, zl, k_, k, i, l, j)
                end do
            end do
        end do
    end do
    ! char(13) is carriage return, so we keep overwriting the percentage
    write (*, "(a1, f5.1,'%')", advance="no") char(13), &
        100._dp * ijkl / size(slater_, 1)
end do
print *
end subroutine

real(dp) function slater_gauss(n, zeta, kk, i, j, k, l) result(r)
! Just like slater(), but uses Gaussian integration to do the integrals.
integer, intent(in) :: n(:)
real(dp), intent(in) :: zeta(:)
integer, intent(in) :: kk, i, j, k, l
real(dp) :: xp
if (n(i)+n(k)-kk-1 < 0 .or. n(j)+n(l)-kk-1 < 0) then
    ! This condition is taken from slater() --- these integrals are
    ! never needed, they don't influence anything, thanks to the wigner3j
    ! selection rule:
    r = 0
else
    r = sto_norm(n(i), zeta(i)) * sto_norm(n(j), zeta(j)) * &
            sto_norm(n(k), zeta(k)) * sto_norm(n(l), zeta(l)) * &
            (gauss_ab(n(i) + n(k), zeta(i) + zeta(k), Ykoverr, 0._dp, 0.1_dp) &
            + gauss_ab(n(i) + n(k), zeta(i) + zeta(k), Ykoverr, 0.1_dp, 5._dp) &
            + spec2(n(i) + n(k), zeta(i) + zeta(k), Ykoverr, 5._dp))
end if

contains

    real(dp) function Ykoverr(x) result(r)
    real(dp), intent(in) :: x
    integer :: n_
    real(dp) :: zeta_
    xp = x ! save x
    n_ = n(j) + n(l)
    zeta_ = zeta(j) + zeta(l)
    r = gauss_ab(n_, zeta_, g, 0._dp, x)
    end function

    real(dp) function g(x)
    real(dp), intent(in) :: x
    !call assert(xp > x)   ! commented out for efficiency reasons
    g = x**kk / xp**(kk+1)
    end function

end function

real(dp) function spec(n, zeta, f) result(res)
! Calculates the integral \int_0^oo r^n * exp(-zeta*r) * f(r) \d r
!
! We convert the integral to:
!
!   \int_0^oo r^n * exp(-zeta*r) * f(r) \d r =
!       = (1/zeta) * \int_0^oo exp(-x) * (x/zeta)^n * f(x/zeta) \d x
!
! And use Gauss-Laguerre quadrature for the integral over "x".
integer, intent(in) :: n
real(dp), intent(in) :: zeta
interface
    real(dp) function f(x)
    import :: dp
    implicit none
    real(dp), intent(in) :: x
    end function
end interface
real(dp) :: r
real(dp), allocatable :: hq(:)
integer :: i

allocate(hq(size(xiq)))
do i = 1, size(xiq)
    r = xiq(i) / zeta
    ! Note: r**n goes to infinity for large "r", and f(r) goes to zero but very
    ! slowly, so r**n * f(r) blows up. For example for Nq=52, the largest point
    ! is x=188.41, zeta=1.7846, so r=105.57, and f(r)=1e-5, but
    ! r**n*f(r)=1837.76
    ! The weight is 1e-80, so the result is zero.
    hq(i) = r**n * f(r)
end do
res = sum(wtq * hq) / zeta
end function

real(dp) function spec2(n, zeta, f, x0) result(res)
! Calculates the integral \int_x0^oo r^n * exp(-zeta*r) * f(r) \d r
!
! We first shift the integral to (0, oo) and then convert the integral to:
!
!   \int_0^oo r^n * exp(-zeta*r) * f(r) \d r =
!       = (1/zeta) * \int_0^oo exp(-x) * (x/zeta)^n * f(x/zeta) \d x
!
! And use Gauss-Laguerre quadrature for the integral over "x".
integer, intent(in) :: n
real(dp), intent(in) :: zeta
interface
    real(dp) function f(x)
    import :: dp
    implicit none
    real(dp), intent(in) :: x
    end function
end interface
real(dp), intent(in) :: x0
real(dp) :: r
real(dp), allocatable :: hq(:)
integer :: i

allocate(hq(size(xiq)))
do i = 1, size(xiq)
    r = xiq(i) / zeta
    ! Note: r**n goes to infinity for large "r", and f(r) goes to zero but very
    ! slowly, so r**n * f(r) blows up. For example for Nq=52, the largest point
    ! is x=188.41, zeta=1.7846, so r=105.57, and f(r)=1e-5, but
    ! r**n*f(r)=1837.76
    ! The weight is 1e-80, so the result is zero.
    hq(i) = (r+x0)**n * f(r+x0)
end do
res = sum(wtq * hq) / zeta
res = res * exp(-zeta*x0)
end function

real(dp) function gauss_ab(n, zeta, f, a, b) result(res)
! Calculates the integral \int_a^b r^n * exp(-zeta*r) * f(r) \d r
!
! A direct Gauss-Legendre quadrature is used.
integer, intent(in) :: n
real(dp), intent(in) :: zeta
interface
    real(dp) function f(x)
    import :: dp
    implicit none
    real(dp), intent(in) :: x
    end function
end interface
real(dp), intent(in) :: a, b
real(dp) :: r
real(dp) :: fq(Nq), jac
integer :: i

jac = (b-a)/2
do i = 1, Nq
    r = (xiq3(i)+1) * jac + a
    fq(i) = r**n * exp(-zeta*r) * f(r)
end do
res = sum(wtq3 * fq * jac)
end function

real(dp) function slater_gauss_screening(n, zeta, kk, i, j, k, l, D) result(r)
! Just like slater(), but uses Gaussian integration to do the integrals.
integer, intent(in) :: n(:)
real(dp), intent(in) :: zeta(:)
integer, intent(in) :: kk, i, j, k, l
real(dp), intent(in) :: D
real(dp) :: xp
if (n(i)+n(k)-kk-1 < 0 .or. n(j)+n(l)-kk-1 < 0) then
    ! This condition is taken from slater() --- these integrals are
    ! never needed, they don't influence anything, thanks to the wigner3j
    ! selection rule:
    r = 0
else
    r = sto_norm(n(i), zeta(i)) * sto_norm(n(j), zeta(j)) * &
            sto_norm(n(k), zeta(k)) * sto_norm(n(l), zeta(l)) * &
            (gauss_ab(n(i) + n(k), zeta(i) + zeta(k), Ykoverr, 0._dp, 0.1_dp) &
            + spec2(n(i) + n(k), zeta(i) + zeta(k), Ykoverr, 0.1_dp))
end if

contains

    real(dp) function Ykoverr(x) result(r)
    real(dp), intent(in) :: x
    integer :: n_
    real(dp) :: zeta_
    xp = x ! save x
    n_ = n(j) + n(l)
    zeta_ = zeta(j) + zeta(l)
    r = gauss_ab(n_, zeta_, g, 0._dp, x)
    end function

    real(dp) function g(x)
    real(dp), intent(in) :: x
    !call assert(xp > x)   ! commented out for efficiency reasons
    g = Vk(kk, D, xp, x)
    end function

end function

subroutine slater_sto_screen(nbfl, nlist, zetalist, slater_, D)
! Just like stoints2, but only the slater integral and uses the Gaussian
! integration to do it. Uses Debye screening.
integer, intent(in) :: nbfl(0:), nlist(:, 0:)
real(dp), intent(in) :: zetalist(:, 0:)
real(dp), allocatable, intent(out) :: slater_(:, :)
real(dp), intent(in) :: D ! Debye screening length
integer,  dimension(sum(nbfl)) :: nl
real(dp), dimension(sum(nbfl)) :: zl
integer :: n, i, j, k, l, k_, ijkl, Lmax, ndof, m
real(dp), allocatable :: dd(:, :)
integer, allocatable :: intindex(:, :, :, :)
! Precalculate factorials. The maximum factorial needed is given by 4*maxn-1:
call calc_factorials(maxval(nlist)*4-1, fact)
print *, "Calculating Slater integrals with screening..."
Lmax = ubound(nbfl, 1)
ndof = sum(nbfl)
m = ndof*(ndof+1)/2
allocate(slater_(m*(m+1)/2, 0:2*Lmax))
slater_ = 0
j = 1
do l = 0, ubound(nbfl, 1)
    nl(j:j+nbfl(l)-1) = nlist   (:nbfl(l), l)
    zl(j:j+nbfl(l)-1) = zetalist(:nbfl(l), l)
    j = j + nbfl(l)
end do

if (.not. allocated(xiq)) then
    !call loadtxt("lag06.txt", d)
    !call loadtxt("lag20.txt", d)
    call loadtxt("lag52.txt", dd)
    allocate(xiq(size(dd, 1)), wtq(size(dd, 1)))
    xiq = dd(:, 1)
    wtq = dd(:, 2)
end if
if (.not. allocated(xiq3)) then
    allocate(xiq3(Nq), wtq3(Nq))
    xiq3 = gauss_pts(Nq)
    wtq3 = gauss_wts(Nq)
end if

n = size(nl)
allocate(intindex(n, n, n, n))
call create_index(intindex)
!$omp parallel default(none) shared(n, nbfl, nl, zl, slater_, D, intindex) private(i, j, k, l, k_, ijkl)
!$omp do schedule(dynamic)
do i = 1, n  ! Only this outer loop is parallelized
    do j = i, n
        do k = 1, n
            do l = k, min((i-1)*i/2 + j - (k-1)*k/2, n)
                ijkl = intindex(l, k, j, i)
                do k_ = 0, 2*ubound(nbfl, 1)
                    !slater_(ijkl, k_) = slater(nl, zl, k_, i, k, j, l)
                    !slater_(ijkl, k_) = slater_gauss(nl, zl, k_, i, k, j, l)
                    !slater_(ijkl, k_) = slater_gauss2(nl, zl, k_, i, k, j, l)
                    !slater_(ijkl, k_) = slater_gauss3(nl, zl, k_, i, k, j, l)
                    slater_(ijkl, k_) = slater_gauss_screening(nl, zl, &
                        k_, i, k, j, l, D) + &
                    slater_gauss_screening(nl, zl, &
                        k_, k, i, l, j, D)
                end do
            end do
        end do
    end do
!    if (omp_get_thread_num() == 0) then
!        ! char(13) is carriage return, so we keep overwriting the percentage
!        write (*, "(a1, f5.1,'%')", advance="no") char(13), 100._dp * i / n
!    end if
end do
!$omp end do
!$omp end parallel
print *
end subroutine

subroutine create_index(intindex)
integer, intent(out) :: intindex(:, :, :, :)
integer :: i, j, k, l, ijkl, n
n = size(intindex, 1)
! Uncomment the commented lines to run checks:
!intindex = -1
ijkl = 1
do i = 1, n
    do j = 1, i
        do k = 1, n
            do l = 1, k
                if ((i-1)*i/2 + j < (k-1)*k/2 + l) cycle
!                call assert(ijkl2intindex(i, j, k, l) == ijkl)
                intindex(i, j, k, l) = ijkl
                intindex(j, i, k, l) = ijkl
                intindex(j, i, l, k) = ijkl
                intindex(i, j, l, k) = ijkl
                intindex(k, l, i, j) = ijkl
                intindex(l, k, i, j) = ijkl
                intindex(l, k, j, i) = ijkl
                intindex(k, l, j, i) = ijkl
                ijkl = ijkl + 1
            end do
        end do
    end do
end do
!call assert(.not. (any(intindex == -1)))
end subroutine

end module
