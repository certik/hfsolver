program ewald_conv4

! Runs a convergence study on the test_ewald4 problem, the example 4.1 in [1].
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
real(dp), parameter :: alpha = 0.90980836237716_dp ! Madelung constant
real(dp), parameter :: stress0 = 1.213077816502878_dp

real(dp) :: stress(6)
real(dp), allocatable :: xred(:, :), forces(:, :), q(:)
real(dp) :: L, E_ewald, E_ewald_fft, E_madelung
integer :: i, Ng, u

natom = 2
ntypat = 1
allocate(xred(3, natom), forces(3, natom), q(natom))
xred(:, 1) = [0, 0, 0]
xred(:, 2) = [1, 1, 1] / 2._dp
q = [1, 1]
L = 2 / sqrt(3._dp)

call ewald_box(L, xred*L, q, E_ewald, forces, stress)
E_ewald = E_ewald / (natom/ntypat)
E_madelung = -2*alpha/L
print *, "a =", L/ang2bohr*100, "pm"
print *, "Madelung: ", E_madelung / kJmol2Ha, "kJ/mol =", E_madelung, "a.u."
print *, "Ewald:    ", E_ewald / kJmol2Ha, "kJ/mol =", E_ewald, "a.u."
print *, "Ewald error:     ", abs(E_ewald - E_madelung), "a.u."
call assert(abs(E_ewald - E_madelung) < 5e-14_dp)
call assert(all(abs(forces) < 1e-15_dp))
call assert(all(abs(stress - (-stress0/L)*[1, 1, 1, 0, 0, 0]) < 1e-15_dp))

open(newunit=u, file="ewald_conv4.txt", status="replace")

Ng = 1
do i = 1, 8
    Ng = Ng*2
    print *
    print *, "Ng =",Ng

    call ewald_fft1(L, xred*L, q, Ng, 0.25_dp*L, E_ewald_fft)
    E_ewald_fft = E_ewald_fft / (natom/ntypat)

    print *, "Ewald FFT:", E_ewald_fft / kJmol2Ha, "kJ/mol =",E_ewald_fft,"a.u."
    print *, "Ewald FFT error: ", abs(E_ewald_fft - E_madelung), "a.u."

    write(u, *) Ng, E_ewald_fft
end do

close(u)

end program
