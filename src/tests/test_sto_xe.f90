program test_sto_xe
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

call get_basis2(54, nbfl, nl, zl, focc)

Z = 54
tolE = 1e-10_dp
tolP = 1e-4_dp
alpha = 0.6_dp
Nscf = 100

n = maxval(nbfl)
Lmax = ubound(nbfl, 1)
ndof = sum(nbfl)
print *, "total  DOFs =", ndof
allocate(S(n, n, 0:Lmax), T(n, n, 0:Lmax), V(n, n, 0:Lmax))
m = ndof*(ndof+1)/2
allocate(slater(m*(m+1)/2, 0:2*Lmax))
call stoints2(Z, nbfl, nl, zl, S, T, V, slater)

allocate(P_(n, n, 0:Lmax), C(n, n, 0:Lmax), lam(n, 0:Lmax))


H = T + V

print *, "SCF cycle:"
call doscf(nbfl, H, slater, S, focc, Nscf, tolE, tolP, alpha, C, P_, lam, Etot)
Ekin = kinetic_energy(nbfl, P_, T)
call printlam(nbfl, lam, Ekin, Etot)

call assert(abs(Etot      - (-7232.1302_dp  )) < 1e-4_dp)
call assert(abs(Ekin      - ( 7232.0468_dp  )) < 1e-4_dp)
call assert(abs(lam(1, 0) - (-1224.397165_dp)) < 1e-6_dp)
call assert(abs(lam(2, 0) - ( -189.340921_dp)) < 1e-6_dp)
call assert(abs(lam(3, 0) - (  -40.175825_dp)) < 1e-6_dp)
call assert(abs(lam(4, 0) - (   -7.856199_dp)) < 1e-6_dp)
call assert(abs(lam(5, 0) - (   -0.944334_dp)) < 1e-6_dp)
call assert(abs(lam(1, 1) - ( -177.782513_dp)) < 1e-6_dp)
call assert(abs(lam(2, 1) - (  -35.221731_dp)) < 1e-6_dp)
call assert(abs(lam(3, 1) - (   -6.008236_dp)) < 1e-6_dp)
call assert(abs(lam(4, 1) - (   -0.457190_dp)) < 1e-6_dp)
call assert(abs(lam(1, 2) - (  -26.119380_dp)) < 1e-6_dp)
call assert(abs(lam(2, 2) - (   -2.777800_dp)) < 1e-6_dp)

end program
