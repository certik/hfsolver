module basis
use types, only: dp
use constants, only: pi
use utils, only: stop_error, lowcase, str, assert
use qc, only: getS, getT, getV, coulomb_repulsion
use gaussians, only: getints2
use scf, only: ijkl2intindex
use basis_aux, only: get_2ints_size
implicit none
private
public normalize, read_basis, get_basis, get_cartesian_shell, getstart, &
        gaussints, get_2ints_size

contains

elemental real(dp) function normalize(l, m, n, alpha) result(norm)
integer, intent(in) :: l, m, n
real(dp), intent(in) :: alpha
norm = sqrt(2**(2*(l+m+n)+1.5)*alpha**(l+m+n+1.5_dp) / &
    fact2(2*l-1)/fact2(2*m-1)/fact2(2*n-1)/pi**(1.5_dp))
end function

pure recursive real(dp) function fact2(n) result(r)
integer, intent(in) :: n
if (n <= 1) then
    r = 1
else
    r = n * fact2(n-2)
end if
end function

subroutine read_basis(Z, lindex, nprimitives, c, zeta)
! Reads the basis from the file
! This function doesn't do any transformation on the basis, so it only returns
! the set of contracted gaussians shells, that are specified on the basis file.
integer, intent(in) :: Z ! atomic number
integer, allocatable, intent(out) :: lindex(:)
! The contracted gaussian L number
integer, allocatable, intent(out) :: nprimitives(:)
! The number of primitive gaussians that compose the given contracted gaussian.
! size(nprimitives) == size(lindex)
real(dp), allocatable, intent(out) :: c(:), zeta(:)
! coefficients and zeta of the primitive gaussians
! size(c) == size(zeta), the first nprimitives(1) correspond to the first
! contracted Gaussian, the next nprimitives(2) correspond to the second
! contracted Gaussian and so on.
integer :: u, ZZ, ios
open(newunit=u, file="p631ss.txt", status="old")
do
    read(u, *, iostat=ios) ZZ
    if (ios /= 0) then
        call stop_error("read_basis: atom Z=" // str(Z) // " not found")
    end if
    call read_atom(u, lindex, nprimitives, c, zeta)
    if (ZZ == Z) then
        close(u)
        return
    end if
end do
end subroutine

subroutine read_atom(u, lindex, nprimitives, c, zeta)
integer, intent(in) :: u
integer, allocatable, intent(out) :: lindex(:), nprimitives(:)
real(dp), allocatable, intent(out) :: c(:), zeta(:)
integer, parameter :: MAX_FCN = 100
real(dp), dimension(MAX_FCN) :: c_tmp, zeta_tmp
character :: L
integer :: r, i, j, idx
read(u, *) r
allocate(lindex(r), nprimitives(r))
idx = 0
do i = 1, r
    read(u, *) L, nprimitives(i)
    lindex(i) = str2l(L)
    do j = 1, nprimitives(i)
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

! Primitive gaussian is described by (l, m, n), norm, coef, alpha,
! (xcenter, ycenter, zcenter).
! Contracted gaussian is described by "nprim" primitive gaussians.
! lambda = l+m+n (0=S, 1=P, 2=D, 3=F, ...) is angular momentum/orbital quantum
! number.
! The basis set file contains pairs of:
! (lambda, contracted gaussian)
! for each atom. For example the 6-31G** basis contains:
!   He: (S, S, P)           # i.e. 1s, 2s, 2p
!   Be: (S, S, P, S, P, D)  # i.e. 1s, 2s, 2p, 3s, 3p, 3d
! To each "S, P, D, F" there is one contracted gaussian, which is described for
! example by:
!           1264.58570000           0.00194480
!            189.93681000           0.01483510
!             43.15908900           0.07209060
!             12.09866300           0.23715420
!              3.80632320           0.46919870
!              1.27289030           0.35652020
! This corresponds to 6 primitive gaussians (alpha, coef) pairs. The
! normalization of *contracted* gaussians has to be computed using the overlap
! matrix (so that it has 1.0 on the diagonal).
! Each contracted gaussian of orbital quantum number "lambda" is actually a
! shell of contracted Gaussians. Two options are used:
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

subroutine get_basis(atZ, at_xyz, nprim, istart, center, power, coef, alpha)
! Returns the basis in a form suitable for calculations.
!
! The input are atomic coordinates and numbers (atZ, at_xyz),
! the output are arrays "nprim, istart, center, power"
! describing each basis function (=contracted gaussian), as well as 'coef' and
! 'alpha' arrays describing the primitive gaussians.
!
integer, intent(in) :: atZ(:) ! atomic numbers
! at_xyz(:, i) are (x,y,z) coordinates of the i-th atom
real(dp), intent(in) :: at_xyz(:, :)
! nprim(i) is the number of primitive gaussians for the i-th contracted
! gaussian
integer, intent(out), allocatable :: nprim(:)
! For the i-th contracted gaussian, the indices i1=istart(i) and
! i2=i1+nprim(i)-1 give all the primitive gaussians: coef(i1:i2) and
! alpha(i1:i2)
integer, intent(out), allocatable :: istart(:)
! center(:, i) are the xyz coordinates of the i-th contracted gaussian
real(dp), intent(out), allocatable :: center(:, :)
! The power(:, i) are the (l,m,n) powers of the i-th contracted gaussian ---
! all functions in the given shell are always given consecutively using the
! get_cartesian_shell() ordering
integer, intent(out), allocatable :: power(:, :)
! The coefficients (including normalization) and alpha's of primitive gaussians
real(dp), intent(out), allocatable :: coef(:), alpha(:)

! Max length of the temporary arrays:
integer, parameter :: N = 10000
integer, allocatable :: at_nprim(:)
integer, allocatable :: at_istart(:)
real(dp), allocatable :: at_center(:, :)
integer, allocatable :: at_power(:, :)
real(dp), allocatable :: at_coef(:), at_alpha(:)
integer :: tmp_nprim(N)
integer :: tmp_istart(N)
real(dp) :: tmp_center(3, N)
integer :: tmp_power(3, N)
real(dp) :: tmp_coef(N), tmp_alpha(N)
integer :: i
integer :: nc ! the number of contracted gaussians (basis functions)
integer :: np ! the number of primitive gaussians
integer :: atc, atp ! number of contracted/primitive gaussians on an atom
nc = 0
np = 0
do i = 1, size(atZ)
    call get_basis_atom(atZ(i), at_xyz(:, i), at_nprim, at_istart, at_center, &
        at_power, at_coef, at_alpha)
    atc = size(at_nprim)
    atp = size(at_coef)
    if (nc+atc > N .or. np+atp > N) then
        call stop_error("get_basis() internal error: Increase N")
    end if
    tmp_nprim(nc+1:nc+atc) = at_nprim
    ! We need to shift the istart array:
    tmp_istart(nc+1:nc+atc) = at_istart + np
    tmp_center(:, nc+1:nc+atc) = at_center
    tmp_power(:, nc+1:nc+atc) = at_power
    tmp_coef(np+1:np+atp) = at_coef
    tmp_alpha(np+1:np+atp) = at_alpha
    nc = nc + atc
    np = np + atp
end do
allocate(nprim(nc), istart(nc), center(3, nc), power(3, nc))
allocate(coef(np), alpha(np))
nprim = tmp_nprim(:nc)
istart = tmp_istart(:nc)
center = tmp_center(:, :nc)
power = tmp_power(:, :nc)
coef = tmp_coef(:np)
alpha = tmp_alpha(:np)
end subroutine

subroutine get_basis_atom(Z, xyz, nprim, istart, center, power, coef, alpha)
! Returns the basis for the atom Z in a form suitable for calculations.
!
! The input are the atom coordinates and number (Z, xyz)
! the output are arrays "nprim, istart, center, power"
! describing each basis function (=contracted gaussian), as well as 'coef' and
! 'alpha' arrays describing the primitive gaussians.
!
integer, intent(in) :: Z ! atomic number
! (x,y,z) coordinates of the i-th atom
real(dp), intent(in) :: xyz(:)
! nprim(i) is the number of primitive gaussians for the i-th contracted
! gaussian
integer, intent(out), allocatable :: nprim(:)
! For the i-th contracted gaussian, the indices i1=istart(i) and
! i2=i1+nprim(i)-1 give all the primitive gaussians: coef(i1:i2) and
! alpha(i1:i2)
integer, intent(out), allocatable :: istart(:)
! center(:, i) are the xyz coordinates of the i-th contracted gaussian
real(dp), intent(out), allocatable :: center(:, :)
! The power(:, i) are the (l,m,n) powers of the i-th contracted gaussian ---
! all functions in the given shell are always given consecutively using the
! get_cartesian_shell() ordering
integer, intent(out), allocatable :: power(:, :)
! The coefficients (including normalization) and alpha's of primitive gaussians
real(dp), intent(out), allocatable :: coef(:), alpha(:)

real(dp), allocatable :: normp(:)
integer, allocatable :: lindex(:), nprimitives(:)
real(dp), allocatable :: c(:), zeta(:)
integer :: n, n2, i, iprim, iprim2, lambda, i1, i2, k

call read_basis(Z, lindex, nprimitives, c, zeta)
n = sum((lindex+1)*(lindex+2)/2)
allocate(nprim(n), istart(n), power(3, n), center(3, n))
do i = 1, n
    center(:, i) = xyz
end do
i = 0
iprim2 = 1
! This do loop creates the arrays describing contracted gaussians (those are
! the actual basis functions) -- there are k = (l+1)*(l+2)/2 contracted
! gaussians in the "Cartesian shell" for the given "l".
do iprim = 1, size(nprimitives)
    i = i + 1
    k = (lindex(iprim)+1)*(lindex(iprim)+2)/2
    nprim(i:i+k-1) = nprimitives(iprim)
    call get_cartesian_shell(lindex(iprim), power(1, i:i+k-1), &
        power(2, i:i+k-1), power(3, i:i+k-1))
    i = i + k - 1
end do
! For each contracted gaussian above, we need to keep track of the individual
! primitive gaussians. Below we create lists of all primitive gaussians and the
! istart() array for accessing it.
istart = getstart(nprim)
n2 = istart(n) + nprim(n) - 1
allocate(normp(n2), coef(n2), alpha(n2))
iprim = 1
i = 1
do
    lambda = sum(power(:, i))
    ! Loop over all "l" for the given Cartesian shell
    do k = 1, (lambda+1)*(lambda+2)/2
        i1 = istart(i)
        i2 = i1 + nprim(i) - 1
        !print *, i, k, iprim
        coef(i1:i2) = c(iprim:iprim+nprim(i)-1)
        alpha(i1:i2) = zeta(iprim:iprim+nprim(i)-1)
        normp(i1:i2) = normalize(power(1, i), power(2, i), power(3, i), &
            alpha(i1:i2))
        i = i + 1
    end do
    iprim = iprim + nprim(i-1)
    if (i > n) exit
end do
coef = coef*normp
call fix_normalization(nprim, istart, center, power, coef, alpha)
end subroutine

subroutine fix_normalization(nprim, istart, center, power, coef, alpha)
! Determines the correct normalization using the overlap integral.
! The result is saved in 'coef'.
integer, intent(in), allocatable :: nprim(:)
integer, intent(in), allocatable :: istart(:)
real(dp), intent(in), allocatable :: center(:, :)
integer, intent(in), allocatable :: power(:, :)
real(dp), intent(inout), allocatable :: coef(:)
real(dp), intent(in), allocatable :: alpha(:)
real(dp) :: S(size(nprim), size(nprim))
integer, dimension(size(power, 2)) :: lpower, mpower, npower
real(dp), dimension(size(center, 2)) :: xcenter, ycenter, zcenter
integer :: i, i1, i2
! The contracted functions must be normalized using the overlap integral,
! currently we calculate the whole S matrix, but we only need the diagonal
! entries.
lpower = power(1, :)
mpower = power(2, :)
npower = power(3, :)
xcenter = center(1, :)
ycenter = center(2, :)
zcenter = center(3, :)
call getS(size(nprim), nprim, istart-1, &
    xcenter, ycenter, zcenter, &
    lpower, mpower, npower, size(coef), coef, alpha, S)
do i = 1, size(nprim)
    if (nprim(i) > 1) then
        ! If contracted, then divide all corresponding primitive gaussians by
        ! the sqrt of the overlap integral
        i1 = istart(i)
        i2 = i1 + nprim(i)-1
        coef(i1:i2) = coef(i1:i2) / sqrt(S(i, i))
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

integer :: i

integer :: n, m
integer, allocatable :: nprim(:)
integer, allocatable :: istart(:)
real(dp), allocatable :: center(:, :)
integer, allocatable :: power(:, :)
real(dp), allocatable :: coef(:), alpha(:)

integer :: Nelec

integer, allocatable :: lpower(:), mpower(:), npower(:)
real(dp), allocatable :: xcenter(:), ycenter(:), zcenter(:)
real(dp), dimension(size(atno)) :: xatom, yatom, zatom

Nelec = sum(atno)

print *, "Creating the basis..."
call get_basis(atno, xyz, nprim, istart, center, power, coef, alpha)
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
