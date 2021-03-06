program test_sto_be5

! The results here are to be compared against [1]:

! [1] Doll, J. D., & Reinhardt, W. P. (1972). Many-Body Green’s Functions for
! Finite, Nonuniform Systems: Applications to Closed Shell Atoms. The Journal
! of Chemical Physics, 57(3), 1169–1184. doi:10.1063/1.1678374

use types, only: dp
use sto, only: stoints2, get_basis2, slater_sto_screen
use utils, only: assert
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

real(dp), allocatable :: S(:, :, :), T(:, :, :), V(:, :, :), slater(:, :)
integer :: n, Z, Nscf, Lmax, ndof
real(dp) :: alpha, Etot, tolE, tolP, Ekin
real(dp), allocatable :: H(:, :, :), P_(:, :, :), C(:, :, :), lam(:, :)
real(dp) :: D

Lmax = 2
allocate(nbfl(0:Lmax), nl(9, 0:Lmax), zl(9, 0:Lmax), focc(2, 0:Lmax))
nbfl = 0
focc = 0
focc(:2, 0) = [2, 2]
nbfl(0) = 5
nl(:5, 0) = [1, 1, 3, 2, 2]
zl(:5, 0) = [5.4297_dp, 2.9954_dp, 3.5810_dp, 1.1977_dp, 0.8923_dp]
nbfl(1) = 5
nl(:5, 1) = [2, 2, 4, 3, 3]
zl(:5, 1) = [5.6998_dp, 2.7850_dp, 4.1500_dp, 1.4387_dp, 0.9819_dp]
nbfl(2) = 2
nl(:2, 2) = [3, 3]
zl(:2, 2) = [1.2662_dp, 7.8314_dp]

Z = 4
D = 20._dp
tolE = 1e-10_dp
tolP = 1e-4_dp
alpha = 0.6_dp
Nscf = 100

n = maxval(nbfl)
ndof = sum(nbfl)
print *, "total  DOFs =", ndof
call stoints2(Z, nbfl, nl, zl, S, T, V, slater)
call slater_sto_screen(nbfl, nl, zl, slater, D)

allocate(P_(n, n, 0:Lmax), C(n, n, 0:Lmax), H(n, n, 0:Lmax), lam(n, 0:Lmax))

H = T + V

print *, "SCF cycle:"
call doscf(nbfl, H, slater, S, focc, Nscf, tolE, tolP, alpha, C, P_, lam, Etot)
Ekin = kinetic_energy(nbfl, P_, T)
call printall(nbfl, nl, zl, lam, C, Ekin, Etot)
print *, "Debye length D =", D
call printlam(nbfl, lam, Ekin, Etot)


call assert(abs(Etot - (-14.85506928_dp)) < 1e-8_dp)
call assert(all(abs(lam(:2, 0) - [-4.86840102_dp, -0.44663327_dp]) < 1e-8_dp))
call assert(all(abs(lam(:2, 1) - [-0.11392548_dp,  0.18282851_dp]) < 1e-8_dp))
call assert(all(abs(lam(:2, 2) - [ 0.40729751_dp, 25.10908153_dp]) < 1e-8_dp))

end program
