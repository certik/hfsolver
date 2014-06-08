program test_ewald2

! This test compares lattice energy 'eew' calculated from the Ewald summation
! against the energy calculated using the Madelung constant for various lattice
! constants L. The agreement is essentially to machine precision.

! Similar to test_ewald, but here we the diagonal Na atom is moved by 3/8
! towards the Cl atom.

use types, only: dp
use constants, only: ang2bohr, kJmol2Ha
use ewald_sums, only: ewald, ewald2, direct_sum, fred2fcart
use utils, only: assert, init_random
implicit none

integer :: natom, ntypat
! Various NaCl lattice constants in A
real(dp), parameter :: Llist(*) = [5.6402_dp, 5.5_dp, 4.5_dp, 6.5_dp, 10._dp]
! The component of the correct force for the given L
real(dp), parameter :: fcorrect(*) = [8.9613425608948444e-2_dp, &
    9.4240310569359220e-2_dp, 0.14077873554188225_dp, &
    6.7473831827766034e-2_dp, 2.8507693947231127e-2_dp]
real(dp), parameter :: stress1_correct(*) = [-0.62770548707869611_dp, &
    -0.64370627058568419_dp, -0.78675210849361421_dp, &
    -0.54467453664942533_dp, -0.35403844882212621_dp]
real(dp), parameter :: stress2_correct(*) = [-0.54226581445122812_dp, &
    -0.55608866303051230_dp, -0.67966392148173749_dp, &
    -0.47053656102581831_dp, -0.30584876466678168_dp]

real(dp) :: ucvol
integer, allocatable :: typat(:)
real(dp) :: gmet(3, 3), rmet(3, 3), rprim(3, 3), gprim(3, 3), stress(6)
real(dp) :: eew
real(dp), allocatable :: xred(:, :), zion(:), grewtn(:, :), fcart(:, :)
real(dp) :: L, alpha, E_ewald, E_madelung
integer :: i, j

! Madelung constant for NaCl, where the diagonal Na atom is moved by 3/8
! towards the Cl atom
alpha = 2.5088837163575413_dp

! Conventional cell:
natom = 8
ntypat = 2
allocate(xred(3, natom), zion(ntypat), grewtn(3, natom), typat(natom), &
    fcart(3, natom))
! Cl^-
xred(:, 1) = [0, 0, 0]
xred(:, 2) = [1, 1, 0] / 2._dp
xred(:, 3) = [1, 0, 1] / 2._dp
xred(:, 4) = [0, 1, 1] / 2._dp
! Na^+
xred(:, 5) = [1, 1, 1] / 2._dp
xred(:, 6) = [1, 0, 0] / 2._dp
xred(:, 7) = [0, 1, 0] / 2._dp
xred(:, 8) = [0, 0, 1] / 2._dp
xred(:, 5:8) = xred(:, 5:8) - spread([3, 3, 3] / 8._dp, 2, 4)
xred(:, 6:8) = xred(:, 6:8) - floor(xred(:, 6:8))
typat = [1, 1, 1, 1, 2, 2, 2, 2]
zion = [-1, +1]

do i = 1, size(Llist)
    L = Llist(i) * ang2bohr
    print *, L

    rmet = 0
    rmet(1, 1) = L**2
    rmet(2, 2) = L**2
    rmet(3, 3) = L**2

    ! gmet = inv(rmet)
    gmet = 0
    gmet(1, 1) = 1/L**2
    gmet(2, 2) = 1/L**2
    gmet(3, 3) = 1/L**2

    ! ucvol = sqrt(det(rmet))
    ucvol = L**3

    ! Reciprocal primitive vectors (without 2*pi) in cartesian coordinates.
    ! gmet = matmul(transpose(gprim), gprim)
    gprim = 0
    gprim(1, 1) = 1 / L
    gprim(2, 2) = 1 / L
    gprim(3, 3) = 1 / L

    ! Real space primitive vectors
    rprim = 0
    rprim(1, 1) = L
    rprim(2, 2) = L
    rprim(3, 3) = L

    call ewald(eew,gmet,grewtn,natom,ntypat,rmet,typat,ucvol,xred,zion)

    E_ewald = eew / (natom/ntypat)
    E_madelung = -2*alpha/L

    call ewald2(gmet,natom,ntypat,rmet,rprim,gprim,stress,typat,ucvol, &
                xred,zion)

    print *, "a =", L/ang2bohr*100, "pm"
    print *, "Madelung:", E_madelung / kJmol2Ha, "kJ/mol"
    print *, "Ewald:   ", E_ewald / kJmol2Ha, "kJ/mol"
    print *, "error:   ", abs(E_ewald - E_madelung), "a.u."
    call assert(abs(E_ewald - E_madelung) < 1e-14_dp)
    stress = -stress * ucvol
    call assert(all(abs(stress - [stress1_correct(i), stress1_correct(i), &
        stress1_correct(i), stress2_correct(i), stress2_correct(i), &
        stress2_correct(i)]) < 1e-10_dp))

    call fred2fcart(fcart, grewtn, gprim)
    do j = 1, 4
        call assert(all(abs(fcart(:, j) - &
            [-fcorrect(i), -fcorrect(i), -fcorrect(i)]) < 1e-10_dp))
    end do
    do j = 5, 8
        call assert(all(abs(fcart(:, j) - &
            [fcorrect(i), fcorrect(i), fcorrect(i)]) < 1e-10_dp))
    end do
end do
deallocate(xred, zion, grewtn, typat, fcart)

print *, "--------"

! Primitive cell (FCC lattice)
natom = 2
ntypat = 2
allocate(xred(3, natom), zion(ntypat), grewtn(3, natom), typat(natom), &
    fcart(3, natom))
! Cl^-
xred(:, 1) = [0, 0, 0]
! Na^+
xred(:, 2) = [1, 1, 1] / 2._dp
xred(:, 2) = xred(:, 2) - [3, 3, 3] / 8._dp
typat = [1, 2]
zion = [-1._dp, 1._dp]

do i = 1, size(Llist)
    L = Llist(i) * ang2bohr

    rmet(1, :) = [2, 1, 1]
    rmet(2, :) = [1, 2, 1]
    rmet(3, :) = [1, 1, 2]
    rmet = rmet * L**2 / 4

    ! gmet = inv(rmet)
    gmet(1, :) = [ 3, -1, -1]
    gmet(2, :) = [-1,  3, -1]
    gmet(3, :) = [-1, -1,  3]
    gmet = gmet / L**2

    ! ucvol = sqrt(det(rmet))
    ucvol = L**3 / 4

    ! Reciprocal primitive vectors (without 2*pi) in cartesian coordinates.
    ! gmet = matmul(transpose(gprim), gprim)
    gprim(:, 1) = [ 1,  1, -1] / L
    gprim(:, 2) = [-1,  1,  1] / L
    gprim(:, 3) = [ 1, -1,  1] / L

    call ewald(eew,gmet,grewtn,natom,ntypat,rmet,typat,ucvol,xred,zion)

    E_ewald = eew / (natom/ntypat)
    E_madelung = -2*alpha/L
    print *, "a =", L/ang2bohr*100, "pm"
    print *, "Madelung:", E_madelung / kJmol2Ha, "kJ/mol"
    print *, "Ewald:   ", E_ewald / kJmol2Ha, "kJ/mol"
    print *, "error:   ", abs(E_ewald - E_madelung), "a.u."
    call assert(abs(E_ewald - E_madelung) < 1e-14_dp)
    call fred2fcart(fcart, grewtn, gprim)
    call assert(all(abs(fcart(:, 1) - &
        [-fcorrect(i), -fcorrect(i), -fcorrect(i)]) < 1e-10_dp))
    call assert(all(abs(fcart(:, 2) - &
        [fcorrect(i), fcorrect(i), fcorrect(i)]) < 1e-10_dp))
end do
deallocate(xred, zion, grewtn, typat, fcart)
end program
