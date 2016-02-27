program test_ewald6

! This test tests direct_sum() for cells with non-zero dipole moment. We use
! ewald_box() as the correct converged results.

use types, only: dp
use constants, only: ang2bohr, kJmol2Ha
use ewald_sums, only: direct_sum, ewald_box, fft_neutralized
use utils, only: assert
implicit none

integer :: natom, ntypat
! Various NaCl lattice constants in A
real(dp), parameter :: Llist(5) = [5.6402_dp, 5.5_dp, 4.5_dp, 6.5_dp, 10._dp]
real(dp) :: stress(6)
real(dp) :: eew
real(dp), allocatable :: xred(:, :), fcart(:, :), q(:), &
    forces(:, :)

real(dp) :: L, alpha, E_ewald, E_madelung, E_direct, rel
integer :: i, j, ncut

alpha = 1.8285774522233_dp ! Madelung constant

! Conventional cell:
natom = 8
ntypat = 2
allocate(xred(3, natom), fcart(3, natom), q(natom), forces(3, natom))
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
q = [-1, -1, -1, -1, 1, 1, 1, 1]*1._dp

! We shift the y-position here, this makes the dipole moment non-zero
xred(2, 2) = 1._dp/4


do i = 1, size(Llist)
    L = Llist(i) * ang2bohr

    call ewald_box(L, xred*L, q, eew, forces, stress)
    E_ewald = eew / (natom/ntypat)

    ncut = 20
    call direct_sum(q, xred*L, L, ncut, eew, fcart)
    E_direct = eew / (natom/ntypat)

    call fft_neutralized(L, xred*L, q, eew, forces, stress)

    E_madelung = -2*alpha/L
    print *, "a =", L/ang2bohr*100, "pm"
    print *, "Madelung:    ", E_madelung / kJmol2Ha, "kJ/mol"
    print *, "Ewald:       ", E_ewald / kJmol2Ha, "kJ/mol"
    print *, "Direct:      ", E_direct / kJmol2Ha, "kJ/mol"
    print *, "Ewald error: ", abs(E_ewald - E_madelung), "a.u."
    print *, "Direct error:", abs(E_direct - E_madelung), "a.u."
    call assert(abs(E_ewald - E_madelung) < 2e-14_dp)
    call assert(abs(E_direct - E_madelung) < 1e-8_dp)
    print *, "Forces errors:"
    do j = 1, natom
        rel = sqrt(sum((fcart(:, j)-forces(:, j))**2)/sum(forces(:, j)**2))
        print *, j, rel
        call assert(rel < 1e-5_dp)
    end do
end do
end program
