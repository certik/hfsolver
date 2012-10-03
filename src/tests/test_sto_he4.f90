program test_sto_he4

! The results here are to be compared against [1]:

! [1] Doll, J. D., & Reinhardt, W. P. (1972). Many-Body Green’s Functions for
! Finite, Nonuniform Systems: Applications to Closed Shell Atoms. The Journal
! of Chemical Physics, 57(3), 1169–1184. doi:10.1063/1.1678374

use types, only: dp
use sto, only: stoints2, get_basis2, slater_sto_gauss
use utils, only: assert
use constants, only: pi, ang2bohr, Ha2eV
use radialscf, only: doscf, kinetic_energy, slater2int22, &
    get_basis, radialC2C, &
    radiallam2lam
use hfprint, only: printall, printlam
use mbpt, only: transform_int2, mbpt2, mbpt3, mbpt4, transform_int22
use gf2, only: find_poles, find_pole_diag, total_energy, plot_poles
use sorting, only: argsort
use scf, only: ijkl2intindex, ijkl2intindex2
implicit none

integer, allocatable :: nl(:, :), nbfl(:)
real(dp), allocatable :: zl(:, :), focc(:, :)

real(dp), allocatable :: S(:, :, :), T(:, :, :), V(:, :, :), slater(:, :), &
    slater2(:, :)
integer :: n, Z, m, Nscf, Lmax, ndof
real(dp) :: alpha, Etot, tolE, tolP, Ekin
real(dp), allocatable :: H(:, :, :), P_(:, :, :), C(:, :, :), lam(:, :)

Lmax = 2
allocate(nbfl(0:Lmax), nl(9, 0:Lmax), zl(9, 0:Lmax), focc(2, 0:Lmax))
nbfl = 0
focc = 0
focc(:1, 0) = [2]
nbfl(0) = 5
nl(:5, 0) = [1, 1, 2, 3, 3]
zl(:5, 0) = [1.4191_dp, 2.5722_dp, 4.2625_dp, 3.9979_dp, 5.4864_dp]
nbfl(1) = 4
nl(:4, 1) = [2, 2, 3, 4]
zl(:4, 1) = [2.5834_dp, 3.6413_dp, 5.5308_dp, 5.7217_dp]
nbfl(2) = 3
nl(:3, 2) = [3, 3, 4]
zl(:3, 2) = [3.6365_dp, 4.8353_dp, 6.9694_dp]

Z = 2
tolE = 1e-10_dp
tolP = 1e-4_dp
alpha = 0.6_dp
Nscf = 100

n = maxval(nbfl)
ndof = sum(nbfl)
print *, "total  DOFs =", ndof
allocate(S(n, n, 0:Lmax), T(n, n, 0:Lmax), V(n, n, 0:Lmax))
m = ndof*(ndof+1)/2
allocate(slater(m*(m+1)/2, 0:2*Lmax))
allocate(slater2(m*(m+1)/2, 0:2*Lmax))
call stoints2(Z, nbfl, nl, zl, S, T, V, slater2)
! We calculate the slater integrals using Gauss-Laguerre quadrature:
call slater_sto_gauss(nbfl, nl, zl, slater)
! This test is crucial: it tests that the numerical and analytical integrals
! agree to very high accuracy:
call assert(all(abs(slater2-slater) < 2e-10_dp))

allocate(P_(n, n, 0:Lmax), C(n, n, 0:Lmax), H(n, n, 0:Lmax), lam(n, 0:Lmax))

H = T + V

print *, "SCF cycle:"
call doscf(nbfl, H, slater, S, focc, Nscf, tolE, tolP, alpha, C, P_, lam, Etot)
Ekin = kinetic_energy(nbfl, P_, T)
!call printall(nbfl, nl, zl, lam, C, Ekin, Etot)
call printlam(nbfl, lam, Ekin, Etot)

call assert(abs(Etot - (-2.86167868_dp)) < 1e-8_dp)
call assert(all(abs(lam(:2, 0) - [-0.91804537_dp, 0.88239287_dp]) < 1e-8_dp))
call assert(all(abs(lam(:2, 1) - [ 1.16674233_dp, 5.07923156_dp]) < 1e-8_dp))
call assert(all(abs(lam(:2, 2) - [ 3.79051952_dp,16.77713044_dp]) < 1e-8_dp))

end program
