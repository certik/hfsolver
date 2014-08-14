program test_sto_ne
use types, only: dp
use sto, only: stoints2, get_basis2
use utils, only: assert
use radialscf, only: doscf, kinetic_energy
use hfprint, only: printlam
implicit none

integer, allocatable :: nl(:, :), nbfl(:)
real(dp), allocatable :: zl(:, :), focc(:, :)

real(dp), allocatable :: S(:, :, :), T(:, :, :), V(:, :, :), slater(:, :)
integer :: n, Z, m, Nscf, Lmax, ndof
real(dp) :: alpha, Etot, tolE, tolP, Ekin
real(dp), allocatable :: H(:, :, :), P_(:, :, :), C(:, :, :), lam(:, :)

call get_basis2(10, nbfl, nl, zl, focc)

Z = 10
tolE = 1e-8_dp
tolP = 1e-4_dp
alpha = 0.6_dp
Nscf = 100

n = maxval(nbfl)
Lmax = ubound(nbfl, 1)
ndof = sum(nbfl)
print *, "total  DOFs =", ndof
allocate(H(n, n, 0:Lmax), S(n, n, 0:Lmax), T(n, n, 0:Lmax), V(n, n, 0:Lmax))
m = ndof*(ndof+1)/2
allocate(slater(m*(m+1)/2, 0:2*Lmax))
call stoints2(Z, nbfl, nl, zl, S, T, V, slater)

allocate(P_(n, n, 0:Lmax), C(n, n, 0:Lmax), lam(n, 0:Lmax))


H = T + V

print *, "SCF cycle:"
call doscf(nbfl, H, slater, S, focc, Nscf, tolE, tolP, alpha, C, P_, lam, Etot)
Ekin = kinetic_energy(nbfl, P_, T)
call printlam(nbfl, lam, Ekin, Etot)

call assert(abs(Etot - (-0.12854705E3_dp)) < 1e-5_dp)
call assert(abs(Ekin - (+0.12854681E3_dp)) < 1e-5_dp)
call assert(abs(lam(1, 0) - (-32.772478_dp)) < 1e-6_dp)
call assert(abs(lam(2, 0) - ( -1.930426_dp)) < 1e-6_dp)
call assert(abs(lam(1, 1) - ( -0.850437_dp)) < 1e-6_dp)

end program
