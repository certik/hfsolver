program test_ewald3

! This test compares lattice energy calculated from the Ewald summation
! against the energy calculated using the Madelung constant for various lattice
! constants L. The agreement is essentially to machine precision.

! This test tests ewald_box().


use types, only: dp
use constants, only: ang2bohr, kJmol2Ha
use ewald_sums, only: ewald_box, fred2fcart, ewald_fft1, min_distance
use utils, only: assert, init_random
implicit none

integer :: natom, ntypat
! Various NaCl lattice constants in A
real(dp), parameter :: Llist(5) = [5.6402_dp, 5.5_dp, 4.5_dp, 6.5_dp, 10._dp]
! The component of the correct force for the given L
! The components of the correct force and stress tensor (multiplied by L)
real(dp), parameter :: fcorrect0 = 10.180280846982683_dp
real(dp), parameter :: stress0 = -4.6601722523551734_dp
real(dp), parameter :: stress1 = -6.6903565769534099_dp
real(dp), parameter :: stress2 = -5.7797035916557675_dp

real(dp) :: stress(6)
real(dp), allocatable :: xred(:, :), forces(:, :), q(:)
real(dp) :: L, alpha, E_ewald, E_ewald_fft, E_madelung
integer :: i, j

alpha = 1.74756459463318219064_dp ! Madelung constant for NaCl

! Conventional cell:
natom = 8
ntypat = 2
allocate(xred(3, natom), forces(3, natom), q(natom))
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
q = [-1, -1, -1, -1, 1, 1, 1, 1]

do i = 1, size(Llist)
    L = Llist(i) * ang2bohr
    call ewald_box(L, xred*L, q, E_ewald, forces, stress)
    E_ewald = E_ewald / (natom/ntypat)

    call ewald_fft1(L, xred*L, q, 64, 0.25_dp*L, E_ewald_fft)
    E_ewald_fft = E_ewald_fft / (natom/ntypat)

    E_madelung = -2*alpha/L

    print *, "a =", L/ang2bohr*100, "pm"
    print *, "Madelung: ", E_madelung / kJmol2Ha, "kJ/mol"
    print *, "Ewald:    ", E_ewald / kJmol2Ha, "kJ/mol"
    print *, "Ewald FFT:", E_ewald_fft / kJmol2Ha, "kJ/mol"
    print *, "error:    ", abs(E_ewald - E_madelung), "a.u."
    print *, "error FFT:", abs(E_ewald_fft - E_madelung), "a.u."
    call assert(abs(E_ewald - E_madelung) < 1e-14_dp)
    call assert(abs(E_ewald_fft - E_madelung) < 1e-10_dp)
    call assert(all(abs(forces) < 1e-17_dp))
    call assert(all(abs(stress - (stress0/L)*[1, 1, 1, 0, 0, 0]) < 1e-15_dp))
end do

print *, "--------"

! Madelung constant for NaCl, where the diagonal Na atom is moved by 3/8
! towards the Cl atom
alpha = 2.5088837163575413_dp

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

print *, "half of xred min distance:", min_distance(xred, 1._dp)/2

do i = 1, size(Llist)
    L = Llist(i) * ang2bohr
    call ewald_box(L, xred*L, q, E_ewald, forces, stress)
    E_ewald = E_ewald / (natom/ntypat)

    call ewald_fft1(L, xred*L, q, 64, 0.10825_dp*L, E_ewald_fft)
    E_ewald_fft = E_ewald_fft / (natom/ntypat)

    E_madelung = -2*alpha/L

    print *, "a =", L/ang2bohr*100, "pm"
    print *, "Madelung: ", E_madelung / kJmol2Ha, "kJ/mol"
    print *, "Ewald:    ", E_ewald / kJmol2Ha, "kJ/mol"
    print *, "Ewald FFT:", E_ewald_fft / kJmol2Ha, "kJ/mol"
    print *, "error:    ", abs(E_ewald - E_madelung), "a.u."
    print *, "error FFT:", abs(E_ewald_fft - E_madelung), "a.u."
    call assert(abs(E_ewald - E_madelung) < 1e-14_dp)
    call assert(abs(E_ewald_fft - E_madelung) < 5e-6_dp)
    call assert(all(abs(stress - [stress1, stress1, stress1, &
        stress2, stress2, stress2]/L) < 1e-10_dp))
    do j = 1, 4
        call assert(all(abs(forces(:, j) - (+fcorrect0/L**2)) < 1e-10_dp))
    end do
    do j = 5, 8
        call assert(all(abs(forces(:, j) - (-fcorrect0/L**2)) < 1e-10_dp))
    end do
end do
end program
