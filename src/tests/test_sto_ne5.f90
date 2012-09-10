program test_sto_ne5

use types, only: dp
use sto, only: stoints2, get_basis2, slater_fe, slater_fe_screen
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
real(dp) :: D

Lmax = 2
allocate(nbfl(0:Lmax), nl(9, 0:Lmax), zl(9, 0:Lmax), focc(2, 0:Lmax))
nbfl = 0
focc = 0
focc(:2, 0) = [2, 2]
focc(:1, 1) = [6]
nbfl(0) = 5
nl(:5, 0) = [1, 1, 2, 2, 2]
zl(:5, 0) = [0.938786D+01, 0.159775D+02, 0.369309D+01, 0.840711D+01, &
    0.240628D+01]
nbfl(1) = 4
nl(:4, 1) = [2, 2, 2, 2]
zl(:4, 1) = [0.244759D+01, 0.449687D+01, 0.166828D+01, 0.934007D+01]
nbfl(2) = 3
nl(:3, 2) = [3, 3, 4]
zl(:3, 2) = [0.395080D+01, 0.120542D+01, 0.745549D+01]

! All zetas are scaled:
!zl(:, 1:) = zl(:, 1:) * 0.952_dp

Z = 10
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
D = 20._dp
call slater_fe_screen(nbfl, nl, zl, slater, D)

allocate(P_(n, n, 0:Lmax), C(n, n, 0:Lmax), H(n, n, 0:Lmax), lam(n, 0:Lmax))

H = T + V

print *, "SCF cycle:"
call doscf(nbfl, H, slater, S, focc, Nscf, tolE, tolP, alpha, C, P_, lam, Etot)
Ekin = kinetic_energy(nbfl, P_, T)
call printall(nbfl, nl, zl, lam, C, Ekin, Etot)
call printlam(nbfl, lam, Ekin, Etot)

call assert(abs(Etot - (-130.72798254_dp)) < 1e-8_dp)
call assert(abs(Ekin - ( 128.62107834_dp)) < 1e-8_dp)
call assert(all(abs(lam(:2, 0) - [-33.19441156_dp, -2.35633119_dp]) < 1e-8_dp))
call assert(all(abs(lam(:2, 1) - [ -1.27607614_dp,  0.41171944_dp]) < 1e-8_dp))
call assert(all(abs(lam(:2, 2) - [  0.13657083_dp,  3.25873422_dp]) < 1e-8_dp))

end program
