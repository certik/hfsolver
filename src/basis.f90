module basis

! Gaussian Type Orbitals (GTO) basis utilities.

! Primitive gaussian is described by (l, m, n), norm, coef, alpha,
! (xcenter, ycenter, zcenter). The normalization is given by norm_prim(), which
! is just a simple formula. The "coef" is used as a coefficient in the
! contracted gaussian:
!
! Contracted gaussian is described by "nprim" primitive gaussians with the same
! (l, m, n). For example with nprim=6, we can have:
!           1264.58570000           0.00194480
!            189.93681000           0.01483510
!             43.15908900           0.07209060
!             12.09866300           0.23715420
!              3.80632320           0.46919870
!              1.27289030           0.35652020
! This corresponds to 6 primitive gaussians (alpha, coef) pairs. Contraction
! can also have just one primitive gaussian (nprim=1).
! Important note: the coefficients above are given relative to *normalized*
! primitive gaussians (this is the standard). So one first has to normalize the
! primitive gaussian, then multiply by the coefficient and only then normalize
! the whole contracted gaussian as described below.
!
! lambda = l+m+n (0=S, 1=P, 2=D, 3=F, ...) is angular momentum/orbital quantum
! number.
! The basis set file contains pairs of (lambda, contracted gaussian)
! for each atom. For example the 6-31G** basis contains:
!   He: (S, S, P)           # i.e. 1s, 2s, 2p
!   Be: (S, S, P, S, P, D)  # i.e. 1s, 2s, 2p, 3s, 3p, 3d
! To each "S, P, D, F" there is one contracted gaussian, which is described for
! example by the (alpha, coef) pairs above.
! Normalization of a contracted gaussian is so that the overlap integral is 1
! (for example the overlap matrix S must have 1.0 on the diagonal).
! More info at [2].
! The contracted gaussian is given by (l, m, n), norm, (xcenter, ycenter,
! zcenter) and "nprim" pairs of (alpha, coef).
!
! Shell is a set of contracted Gaussians for all possible combinations of (l, m
! n) for the given lambda. All other parameters stay the same, so the shell is
! given by: lambda, norm, (xcenter, ycenter, zcenter) and "nprim" pairs of
! (alpha, coef). There are two possibilities:
!
! 1) Cartesian Gaussians shell: there are (lambda+1)*(lambda+2)/2 contracted
! gaussians (all combinations of (n, l, m) so that lambda=n+l+m)
!
! 2) Spherical Gaussians shell: there are 2*lambda+1 contracted gaussians (all
! combinations of "m" for the given l=lambda for the spherical harmonic Y_lm)
!
! l    2l+1     (l+1)(l+1)/2   =decomposition in l      sum of 2l+1
! 0s     1           1              0                   1            =  1
! 1p     3           3              1                     3          =  3
! 2d     5           6              0,2                 1 + 5        =  6
! 3f     7          10              1,3                   3 + 7      = 10
! 4g     9          15              0,2,4               1 + 5 + 9    = 15
! 5h    11          21              1,3,5                 3 + 7 + 11 = 21
!
! So for l = 0, 1 the number of basis functions is the same, but for l >= 2 the
! Cartesian Gaussians have more functions. The homogeneous polynomials of order
! "n" in (x, y, z) span the symmetric (n+1)(n+2)/2-dimensional reducible
! representation of the U(3) group. The spherical harmonics Y_lm span the 2l+1
! dimensional irreducible representation of the SO(3) group, which is a
! subgroup of U(3).  From the table above, s, p is the same, cartesian d is
! composed of five spherical d and one s (total of 6 functions), cartesian h is
! composed of three p, seven f and eleven h (total 21 functions).  See [1] for
! more info.
!
! [1] Schlegel, H. B., & Frisch, M. J. (1995). Transformation Between Cartesian
! and Pure Spherical Harmonic Gaussians. International Journal of Quantum
! Chemistry, 54(2), 83-87.
!
! [2] http://www.files.chem.vt.edu/chem-dept/valeev/docs/ints.pdf


! In this file, we only use Cartesian GTO and each shell has exactly
! (lambda+1)*(lambda+2)/2 contracted GTO. The ordering of the (n, l, m)
! triplets is given by the get_cartesian_shell() function.


use types, only: dp
use constants, only: pi
use utils, only: stop_error, lowcase, str, assert
use qc, only: getS, getT, getV, coulomb_repulsion
use gaussians, only: getints2
use scf, only: ijkl2intindex
use basis_aux, only: get_2ints_size
implicit none
private
public get_cartesian_shell, getstart, gaussints, get_2ints_size

contains

elemental real(dp) function norm_prim(l, m, n, alpha) result(norm)
integer, intent(in) :: l, m, n
real(dp), intent(in) :: alpha
norm = sqrt(2**(2*(l+m+n)+1.5)*alpha**(l+m+n+1.5_dp) / &
    fact2(2*l-1)/fact2(2*m-1)/fact2(2*n-1)/pi**(1.5_dp))
end function

pure real(dp) function norm_contr(l, m, n, alpha, coefu) result(norm)
! Calculates the normalization coefficient of a contracted GTO specified by (l,
! m, n) and (alpha, coefu) pairs. The "coefu" here is given with regards to
! *unnormalized* primitive GTOs (thus the "u" in "coefu"), so the normalized
! contracted GTO is then given by:
!     norm_contr * sum_i coefu(i) * GTO(l, m, n, alpha(i))
!
! Note: if there is only one primitive GTO in the contraction (that is,
! size(alpha)==size(coefu)==1), then norm_contr() is equivalent to
! norm_prim(), that is, norm_contr = norm_prim / coefu(1):
!
!    norm_contr * coefu(1) * GTO(l, m, n, alpha(1)) =
!       = norm_prim / coefu(1) * coefu(1) * GTO(l, m, n, alpha(1)) =
!       = norm_prim * GTO(l, m, n, alpha(1))
!
! (l, m, n) of the contracted GTO:
integer, intent(in) :: l, m, n
! Alpha exponents of primitive GTOs:
real(dp), intent(in) :: alpha(:)
! Coefficients with regards to *unnormalized* primitive GTOs:
real(dp), intent(in) :: coefu(:)
real(dp) :: L_, S
integer :: i, j
L_ = l + m + n
S = 0
do i = 1, size(alpha)
    do j = 1, size(alpha)
        S = S + coefu(i)*coefu(j)/(alpha(i)+alpha(j))**(L_+3._dp/2)
    end do
end do
norm = 1/pi**(3._dp/4) * sqrt(2**L_/ fact2(2*l-1)/fact2(2*m-1)/fact2(2*n-1)) / &
    sqrt(S)
end function

pure recursive real(dp) function fact2(n) result(r)
integer, intent(in) :: n
if (n <= 1) then
    r = 1
else
    r = n * fact2(n-2)
end if
end function

subroutine read_basis(Z, L, nprim, c, zeta)
! Reads the basis corresponding to the atom Z from the file.
! Currently it reads the 6-31G** basis.
! This function doesn't do any transformation on the basis, so it only returns
! the set of contracted GTO shells, that are specified in the basis file.

integer, intent(in) :: Z ! atomic number

! The contracted gaussian L number. L(i) is the L number of the i-th shell:
integer, allocatable, intent(out) :: L(:)

! nprim(i) is the number of primitive GTO that compose the given contracted GTO
! of the i-th shell. size(nprim) == size(L)
integer, allocatable, intent(out) :: nprim(:)

! Coefficients and zeta of the primitive gaussians.
! size(c) == size(zeta), the first nprim(1) correspond to the first
! contracted Gaussian, the next nprim(2) correspond to the second
! contracted Gaussian and so on.
real(dp), allocatable, intent(out) :: c(:), zeta(:)
!
! Important note: the coefficients c(:) are with respect to *normalized*
! primitive GTO.
!
! Example:
! --------
!
! Given the part of input file for Lithium (Z=3):
!
! 3
! 6
! S 6
!             642.41892000           0.00214260
!              96.79851500           0.01620890
!              22.09112100           0.07731560
!               6.20107030           0.24578600
!               1.93511770           0.47018900
!               0.63673580           0.34547080
! S 3
!               2.32491840          -0.03509170
!               0.63243060          -0.19123280
!               0.07905340           1.08398780
! P 3
!               2.32491840           0.00894150
!               0.63243060           0.14100950
!               0.07905340           0.94536370
! S 1
!               0.03596200           1.00000000
! P 1
!               0.03596200           1.00000000
! D 1
!               0.20000000           1.00000000
!
! There are 6 shells and the output arrays will be:
!
! L     = [0, 0, 1, 0, 1, 2]
! nprim = [6, 3, 3, 1, 1, 1]
! zeta  = [642.41892000, 96.79851500, ..., 0.03596200, 0.20000000]
! c     = [0.00214260, 0.01620890, ..., 1.00000000, 1.00000000]
!
! The coefficients c(:) are with respect to normalized GTOs. From this it
! follows that primitive GTO corresponding to nprim(i)=1 will have a
! coefficient equal to 1.0. In order to obtain coefficients with respect to
! unnormalized GTOs, one has to multiply each c(:) by the norm factor (in
! general there will be no "1.0"s).

integer :: u, ZZ, ios
open(newunit=u, file="p631ss.txt", status="old") ! 6-31G**
do
    read(u, *, iostat=ios) ZZ
    if (ios /= 0) then
        call stop_error("read_basis: atom Z=" // str(Z) // " not found")
    end if
    call read_atom(u, L, nprim, c, zeta)
    if (ZZ == Z) then
        close(u)
        return
    end if
end do
end subroutine

subroutine read_atom(u, L, nprim, c, zeta)
! Reads the basis of one atom from the file "u".
! The meaning of L, nprim, c and zeta is the same as in read_basis().
integer, intent(in) :: u
integer, allocatable, intent(out) :: L(:), nprim(:)
real(dp), allocatable, intent(out) :: c(:), zeta(:)
integer, parameter :: MAX_FCN = 100
real(dp), dimension(MAX_FCN) :: c_tmp, zeta_tmp
character :: Ltmp
integer :: r, i, j, idx
read(u, *) r
allocate(L(r), nprim(r))
idx = 0
do i = 1, r
    read(u, *) Ltmp, nprim(i)
    L(i) = str2l(Ltmp)
    do j = 1, nprim(i)
        idx = idx + 1
        if (idx > MAX_FCN) call stop_error("read_atom: increase MAX_FCN")
        read(u, *) zeta_tmp(idx), c_tmp(idx)
    end do
end do
allocate(c(idx), zeta(idx))
c = c_tmp(:idx)
zeta = zeta_tmp(:idx)
end subroutine

integer function str2l(s) result(l)
character, intent(in) :: s
select case(lowcase(s))
    case ("s")
        l = 0
    case ("p")
        l = 1
    case ("d")
        l = 2
    case ("f")
        l = 3
    case ("g")
        l = 4
    case ("h")
        l = 5
    case ("i")
        l = 6
    case default
        l = -1 ! Suppress compiler warning
        call stop_error("str2l: unknown angular momentum '" // s // "'")
end select
end function


subroutine get_basis(atZ, at_xyz, L, nprim, center, coef, alpha)
! Returns the basis as a list of *shells*.

! The input are atomic coordinates and numbers (atZ, at_xyz),
! The output are arrays "L(i), nprim(i), center(i)"
! describing i-th shell. The arrays "coef, alpha" contain the primitive GTOs of
! all shells together in consecutive order. In order to obtain coef, alpha
! for the i-th shell, do:
!     istart = getstart(nprim)
!     print *,  coef(istart(i):istart(i)+nprim(i)-1)
!     print *, alpha(istart(i):istart(i)+nprim(i)-1)

integer, intent(in) :: atZ(:) ! atZ(i) is the atomic number of the i-th atom
! at_xyz(:, i) are (x,y,z) coordinates of the i-th atom
real(dp), intent(in) :: at_xyz(:, :)
! L(i) is the lambda of the i-th shell
integer, intent(out), allocatable :: L(:)
! nprim(i) is the number of primitive GTO in the i-th shell:
integer, intent(out), allocatable :: nprim(:)
! center(:, i) are the xyz coordinates of the i-th shell
real(dp), intent(out), allocatable :: center(:, :)
! The coefficients (with regards to normalized GTO) and alpha's of primitive
! gaussians:
real(dp), intent(out), allocatable :: coef(:), alpha(:)

! Max length of the temporary arrays:
integer, parameter :: N = 10000
integer, allocatable :: at_nprim(:), at_L(:)
real(dp), allocatable :: at_center(:, :)
real(dp), allocatable :: at_coef(:), at_alpha(:)
integer :: tmp_L(N)
integer :: tmp_nprim(N)
real(dp) :: tmp_center(3, N)
real(dp) :: tmp_coef(N), tmp_alpha(N)
integer :: i, j
integer :: nc, np, atc, atp
nc = 0
np = 0
do i = 1, size(atZ)
    call read_basis(atZ(i), at_L, at_nprim, at_coef, at_alpha)
    allocate(at_center(3, size(at_L)))
    do j = 1, size(at_L)
        at_center(:, j) = at_xyz(:, i)
    end do

    atc = size(at_nprim)
    atp = size(at_coef)
    if (nc+atc > N .or. np+atp > N) then
        call stop_error("get_basis() internal error: Increase N")
    end if
    tmp_L(nc+1:nc+atc) = at_L
    tmp_nprim(nc+1:nc+atc) = at_nprim
    tmp_center(:, nc+1:nc+atc) = at_center

    tmp_coef(np+1:np+atp) = at_coef
    tmp_alpha(np+1:np+atp) = at_alpha

    nc = nc + atc
    np = np + atp
    deallocate(at_center)
end do
allocate(L(nc), nprim(nc), center(3, nc))
allocate(coef(np), alpha(np))
L = tmp_L(:nc)
nprim = tmp_nprim(:nc)
center = tmp_center(:, :nc)

coef = tmp_coef(:np)
alpha = tmp_alpha(:np)
end subroutine

subroutine fix_normalization(nprim, istart, power, coefu, alpha)
! Determines the correct normalization using the overlap integral.
! The result is saved in 'coefu'.
integer, intent(in), allocatable :: nprim(:)
integer, intent(in), allocatable :: istart(:)
integer, intent(in), allocatable :: power(:, :)
! Coefficients with regards to *unnormalized* primitive gaussians
real(dp), intent(inout), allocatable :: coefu(:)
real(dp), intent(in), allocatable :: alpha(:)
integer, dimension(size(power, 2)) :: l, m, n
integer :: i, i1, i2
real(dp) :: norm
l = power(1, :)
m = power(2, :)
n = power(3, :)
do i = 1, size(nprim)
    if (nprim(i) > 1) then
        i1 = istart(i)
        i2 = i1 + nprim(i)-1
        norm = norm_contr(l(i), m(i), n(i), alpha(i1:i2), coefu(i1:i2))
        coefu(i1:i2) = coefu(i1:i2) * norm
    end if
end do
end subroutine

subroutine get_cartesian_shell(L, lpower, mpower, npower)
! Returns all the functions on the "L" Cartesian shell.
! The Cartesian shell is defined as (x^l, y^m, z^n) where the powers (l, m, n)
! are given by all combinations of l,m,n <= L so that l+m+n=L. There
! is exactly (L+1)(L+2)/2 such combinations.
! There are several ways how one can order these combinations. This subroutine
! uses the following ordering:
! L = 0
! lpower == [0]
! mpower == [0]
! npower == [0]
! L = 1
! lpower == [1, 0, 0]
! mpower == [0, 1, 0]
! npower == [0, 0, 1]
! L = 2
! lpower == [2, 1, 1, 0, 0, 0]
! mpower == [0, 1, 0, 2, 1, 0]
! npower == [0, 0, 1, 0, 1, 2]
! L = 3
! lpower == [3, 2, 2, 1, 1, 1, 0, 0, 0, 0]
! mpower == [0, 1, 0, 2, 1, 0, 3, 2, 1, 0]
! npower == [0, 0, 1, 0, 1, 2, 0, 1, 2, 3]
! ...

integer, intent(in) :: L  ! The Gaussian orbital number
integer, intent(out), dimension(:) :: lpower, mpower, npower
! the (l, m, n) powers. The arrays must have dimension (L+1)(L+2)/2
integer :: i, k, n
i = 1
do k = 0, L
    do n = 0, k
        lpower(i) = L - k
        mpower(i) = k - n
        npower(i) = n
        i = i + 1
    end do
end do
end subroutine

function getstart(nprim) result(istart)
! Calculates the start of the primitive gaussians for the given shell ("nprim")
! Example:
! nprim  = [6, 3, 3, 1, 1, 1]
! istart = [1, 7, 10, 13, 14, 15]
integer, intent(in) :: nprim(:)
integer :: istart(size(nprim))
integer :: i
istart(1) = 1
do i = 2, size(nprim)
    istart(i) = istart(i-1) + nprim(i-1)
end do
end function

subroutine gaussints(atno, xyz, S, T, V, int2)
integer, intent(in) :: atno(:) ! Atomic numbers (Z)
real(dp), intent(in) :: xyz(:, :) ! Atomic coordinates (3, :)
real(dp), allocatable, intent(out) :: S(:, :), & ! Overlap matrix
    T(:, :), &  ! Kinetic matrix
    V(:, :), &  ! Potential matrix
    int2(:)     ! Two particle integrals indexed using ijkl2intindex()

integer :: Nelec

integer, allocatable :: L(:)
integer, allocatable :: nprim(:)
integer, allocatable :: istart(:)
real(dp), allocatable :: coef(:), alpha(:)

integer :: i
integer :: n, m

real(dp), allocatable :: center(:, :)
integer, allocatable :: power(:, :)


integer, allocatable :: lpower(:), mpower(:), npower(:)
real(dp), allocatable :: xcenter(:), ycenter(:), zcenter(:)
real(dp), dimension(size(atno)) :: xatom, yatom, zatom

Nelec = sum(atno)

print *, "Creating the basis..."
call get_basis(atno, xyz, L, nprim, center, coef, alpha)
n = size(nprim)
print "(a)", "----------------------------------------------------------"
print "(a)", "Center     Atomic              Coordinates (a.u.)"
print "(a)", "Number     Number             X           Y           Z"
print "(a)", "----------------------------------------------------------"
do i = 1, size(atno)
    print "(i4, i11, '       ', 3f12.6)", i, atno(i), xyz(:, i)
end do
print "(a)", "----------------------------------------------------------"
print *, "Nelec =", Nelec
print *, "DOFs  =", n
if (n > 120) call stop_error("Would need too much memory.")
m = n*(n+1)/2
allocate(S(n, n), T(n, n), V(n, n))
allocate(int2(m*(m+1)/2))
allocate(lpower(n), mpower(n), npower(n), xcenter(n), ycenter(n), zcenter(n))
lpower = power(1, :)
mpower = power(2, :)
npower = power(3, :)
xcenter = center(1, :)
ycenter = center(2, :)
zcenter = center(3, :)
xatom = xyz(1, :)
yatom = xyz(2, :)
zatom = xyz(3, :)


print *, "Calculating S ..."
call getS(n, nprim, istart-1, xcenter, ycenter, zcenter, &
        lpower, mpower, npower, size(coef), coef, alpha, S)

print *, "Calculating T ..."
call getT(n, nprim, istart-1, xcenter, ycenter, zcenter, &
        lpower, mpower, npower, size(coef), coef, alpha, T)

print *, "Calculating V ..."
call getV(n, nprim, istart-1, xcenter, ycenter, zcenter, &
        lpower, mpower, npower, size(coef), coef, alpha, &
        size(atno), atno, xatom, yatom, zatom, V)

print *, "Calculating two particle integrals (", size(int2), ")..."
call getints2(nprim, istart, center, power, coef, alpha, int2)

end subroutine

end module
