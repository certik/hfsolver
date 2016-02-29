program test_ewald6

! This test tests direct_sum() for cells with non-zero dipole moment. We use
! ewald_box() as the correct converged results.

use types, only: dp
use ewald_sums, only: direct_sum, ewald_box, fft_neutralized
use utils, only: assert
implicit none

integer :: natom, ntypat
! Various NaCl lattice constants in A
real(dp), parameter :: Llist(1) = [2/sqrt(3._dp)]
real(dp) :: stress(6)
real(dp) :: eew
real(dp), allocatable :: xred(:, :), fcart(:, :), q(:), &
    forces(:, :)

real(dp) :: L, alpha, E_ewald, E_madelung, E_direct, E_fft
integer :: i, ncut

alpha = 1.8285774522233_dp ! Madelung constant

! Conventional cell:
natom = 3
ntypat = 2
allocate(xred(3, natom), fcart(3, natom), q(natom), forces(3, natom))
! Na^+
xred(:, 1) = [0._dp, 0._dp, 0._dp]
xred(:, 2) = [1._dp/2, 1._dp/2, 1._dp/2]
xred(:, 3) = [1._dp/4, 1._dp/2, 1._dp/2]
q = [2, -1, -1]*1._dp
!q = [1, 1]*1._dp


do i = 1, size(Llist)
    L = Llist(i)

    call ewald_box(L, xred*L, q, eew, forces, stress)
    E_ewald = eew / (natom/ntypat)

    ncut = 20
    call direct_sum(q, xred*L, L, ncut, eew, fcart)
    E_direct = eew / (natom/ntypat)

    call fft_neutralized(L, xred*L, q, eew)
    E_fft = eew / (natom/ntypat)

    E_madelung = -2*alpha/L
    print *, "a =", L, "a.u."
    print *, "Ewald:       ", E_ewald, "Hartree/atom"
    print *, "FFT:         ", E_fft, "Hartree/atom"
    print *, "Ewald error:", abs(E_ewald - (-1.5758343085)), "a.u."
    print *, "FFT error:", abs(E_ewald - E_fft), "a.u."
end do
end program
