program test_ewald4

! This test compares lattice energy calculated from the Ewald summation
! against the energy calculated using the Madelung constant for various lattice
! constants L. The agreement is essentially to machine precision.

! This test is testing the example from section 4.1 in [1].
!
! [1] Pask, J. E., Sukumar, N., & Mousavi, S. E. (2012). Linear scaling
! solution of the all-electron Coulomb problem in solids. International Journal
! for Multiscale Computational Engineering, 10(1), 83–99.
! doi:10.1615/IntJMultCompEng.2011002201


use types, only: dp
use constants, only: ang2bohr, kJmol2Ha
use ewald_sums, only: ewald_box, ewald_fft1, fred2fcart
use utils, only: assert, init_random
implicit none

integer :: natom, ntypat
real(dp), parameter :: Llist(4) = [2 / sqrt(3._dp), 1._dp, 2._dp, 3._dp]
real(dp), parameter :: alpha = 0.90980836237716_dp ! Madelung constant
real(dp), parameter :: stress0 = 1.213077816502878_dp

real(dp) :: stress(6)
real(dp), allocatable :: xred(:, :), forces(:, :), q(:)
real(dp) :: L, E_ewald, E_ewald_fft, E_madelung
integer :: i

natom = 2
ntypat = 1
allocate(xred(3, natom), forces(3, natom), q(natom))
xred(:, 1) = [0, 0, 0]
xred(:, 2) = [1, 1, 1] / 2._dp
q = [1, 1]

do i = 1, size(Llist)
    L = Llist(i)
    call ewald_box(L, xred*L, q, E_ewald, forces, stress)
    E_ewald = E_ewald / (natom/ntypat)

    call ewald_fft1(L, xred*L, q, 64, 0.25_dp*L, E_ewald_fft)
    E_ewald_fft = E_ewald_fft / (natom/ntypat)

    E_madelung = -2*alpha/L

    print *, "a =", L/ang2bohr*100, "pm"
    print *, "Madelung: ", E_madelung / kJmol2Ha, "kJ/mol =", E_madelung, "a.u."
    print *, "Ewald:    ", E_ewald / kJmol2Ha, "kJ/mol =", E_ewald, "a.u."
    print *, "Ewald FFT:", E_ewald_fft / kJmol2Ha, "kJ/mol =",E_ewald_fft,"a.u."
    print *, "Ewald error:     ", abs(E_ewald - E_madelung), "a.u."
    print *, "Ewald FFT error: ", abs(E_ewald_fft - E_madelung), "a.u."
    call assert(abs(E_ewald - E_madelung) < 5e-14_dp)
    call assert(abs(E_ewald_fft - E_madelung) < 1e-10_dp)
    call assert(all(abs(forces) < 1e-15_dp))
    call assert(all(abs(stress - (-stress0/L)*[1, 1, 1, 0, 0, 0]) < 1e-15_dp))
end do

end program
