program test_ewald

! This test compares lattice energy 'eew' calculated from the Ewald summation
! against the energy calculated using the Madelung constant for various lattice
! constants L. The agreement is essentially to machine precision.

use types, only: dp
use constants, only: ang2bohr, kJmol2Ha
use ewald_sums, only: ewald, fred2fcart, direct_sum
use utils, only: assert, init_random
implicit none

integer :: natom, ntypat
! Various NaCl lattice constants in A
real(dp), parameter :: Llist(*) = [5.6402_dp, 5.5_dp, 4.5_dp, 6.5_dp, 10._dp]
real(dp) :: ucvol
integer, allocatable :: typat(:)
real(dp) :: gmet(3, 3), rmet(3, 3), gprim(3, 3)
real(dp) :: eew
real(dp), allocatable :: xred(:, :), zion(:), grewtn(:, :), fcart(:, :)

real(dp) :: L, alpha, E_ewald, E_madelung, E_direct
integer :: i, ncut

alpha = 1.74756459463318219064_dp ! Madelung constant for NaCl

! Conventional cell:
natom = 8
ntypat = 2
allocate(xred(3, natom), zion(ntypat), grewtn(3, natom), typat(natom), &
    fcart(3, natom))
! Cl^-
xred(:, 1) = [0._dp, 0._dp, 0._dp]
xred(:, 2) = [1._dp/2, 1._dp/2, 0._dp]
xred(:, 3) = [1._dp/2, 0._dp, 1._dp/2]
xred(:, 4) = [0._dp, 1._dp/2, 1._dp/2]
! Na^+
xred(:, 5) = [1._dp/2, 1._dp/2, 1._dp/2]
xred(:, 6) = [1._dp/2, 0._dp, 0._dp]
xred(:, 7) = [0._dp, 1._dp/2, 0._dp]
xred(:, 8) = [0._dp, 0._dp, 1._dp/2]
typat = [1, 1, 1, 1, 2, 2, 2, 2]
zion = [-1._dp, +1._dp]

do i = 1, size(Llist)
    L = Llist(i) * ang2bohr

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

    call ewald(eew,gmet,grewtn,natom,ntypat,rmet,typat,ucvol,xred,zion)

    E_ewald = eew / (natom/ntypat)

    ncut = 20
    call direct_sum([-1, -1, -1, -1, 1, 1, 1, 1]*1._dp, xred*L, L, ncut, eew, &
        fcart)
    E_direct = eew / (natom/ntypat)

    E_madelung = -2*alpha/L
    print *, "a =", L/ang2bohr*100, "pm"
    print *, "Madelung:    ", E_madelung / kJmol2Ha, "kJ/mol"
    print *, "Ewald:       ", E_ewald / kJmol2Ha, "kJ/mol"
    print *, "Direct:      ", E_direct / kJmol2Ha, "kJ/mol"
    print *, "Ewald error: ", abs(E_ewald - E_madelung), "a.u."
    print *, "Direct error:", abs(E_direct - E_madelung), "a.u."
    call fred2fcart(fcart, grewtn, gprim)
    call assert(abs(E_ewald - E_madelung) < 1e-14_dp)
    call assert(abs(E_direct - E_madelung) < 1e-8_dp)
    call assert(all(abs(fcart) < 1e-15_dp))
end do
deallocate(xred, zion, grewtn, typat, fcart)

print *, "--------"

! Primitive cell (FCC lattice)
natom = 2
ntypat = 2
allocate(xred(3, natom), zion(ntypat), grewtn(3, natom), typat(natom), &
    fcart(3, natom))
! Cl^-
xred(:, 1) = [0._dp, 0._dp, 0._dp]
! Na^+
xred(:, 2) = [1._dp/2, 1._dp/2, 1._dp/2]
typat = [1, 2]
zion = [-1._dp, +1._dp]

do i = 1, size(Llist)
    L = Llist(i) * ang2bohr

    rmet(1, :) = [2._dp, 1._dp, 1._dp]
    rmet(2, :) = [1._dp, 2._dp, 1._dp]
    rmet(3, :) = [1._dp, 1._dp, 2._dp]
    rmet = rmet * L**2 / 4

    ! gmet = inv(rmet)
    gmet(1, :) = [3._dp, -1._dp, -1._dp]
    gmet(2, :) = [-1._dp, 3._dp, -1._dp]
    gmet(3, :) = [-1._dp, -1._dp, 3._dp]
    gmet = gmet / L**2

    ! ucvol = sqrt(det(rmet))
    ucvol = L**3 / 4

    ! Reciprocal primitive vectors (without 2*pi) in cartesian coordinates.
    ! gmet = matmul(transpose(gprim), gprim)
    gprim(:, 1) = [1, 1, -1] / L
    gprim(:, 2) = [-1, 1, 1] / L
    gprim(:, 3) = [1, -1, 1] / L

    call ewald(eew,gmet,grewtn,natom,ntypat,rmet,typat,ucvol,xred,zion)

    E_ewald = eew / (natom/ntypat)
    E_madelung = -2*alpha/L
    print *, "a =", L/ang2bohr*100, "pm"
    print *, "Ewald:   ", E_ewald / kJmol2Ha, "kJ/mol"
    print *, "Madelung:", E_madelung / kJmol2Ha, "kJ/mol"
    print *, "error:   ", abs(E_ewald - E_madelung), "a.u."
    call assert(abs(E_ewald - E_madelung) < 1e-14_dp)

    call fred2fcart(fcart, grewtn, gprim)
    call assert(all(abs(fcart) < 1e-15_dp))
end do
deallocate(xred, zion, grewtn, typat, fcart)
end program
