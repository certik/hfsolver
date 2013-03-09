program test_sto_mg
use types, only: dp
use sto, only: stoints2, get_basis2
use utils, only: assert
use radialscf, only: doscf, kinetic_energy
use hfprint, only: printall, printlam
implicit none

integer, allocatable :: nl(:, :), nbfl(:)
real(dp), allocatable :: zl(:, :), focc(:, :)

real(dp), allocatable :: S(:, :, :), T(:, :, :), V(:, :, :), slater(:, :)
integer :: n, Z, m, Nscf, Lmax, ndof
real(dp) :: alpha, Etot, tolE, tolP, Ekin
real(dp), allocatable :: H(:, :, :), P_(:, :, :), C(:, :, :), lam(:, :)

allocate(nbfl(0:2), nl(8, 0:2), zl(8, 0:2), focc(3, 0:2))
nbfl = 0
focc = 0
focc(:3, 0) = [2, 2, 2]
focc(:1, 1) = [6]
nbfl(0) = 8
nl(:8, 0) = [1, 2, 2, 3, 3, 3, 3, 4]
zl(:8, 0) = [12._dp, 13.6_dp, 9.3_dp, 6.5_dp, 4.2_dp, 1.4_dp, 0.9_dp, 2.5_dp]
nbfl(1) = 8
nl(:8, 1) = [2, 2, 2, 3, 3, 3, 3, 4]
zl(:8, 1) = [12.5_dp, 9.2_dp, 6._dp, 0.5_dp, 5.3_dp, 3.7_dp, 2.5_dp, 1.2_dp]
nbfl(2) = 8
nl(:8, 2) = [3, 3, 3, 3, 4, 4, 4, 4]
zl(:8, 2) = [12.8_dp, 9.6_dp, 6.4_dp, 0.75_dp, 5.2_dp, 3.6_dp, 2.1_dp, 1.4_dp]

Z = 12
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
call printall(nbfl, nl, zl, lam, C, Ekin, Etot)
call printlam(nbfl, lam, Ekin, Etot)

call assert(abs(Etot - (-199.61259393_dp)) < 1e-8_dp)
call assert(abs(Ekin - (+199.62293919_dp)) < 1e-8_dp)
call assert(all(abs(lam(:3, 0) - [-49.029177_dp, -3.765885_dp, &
    -0.252718_dp]) < 1e-6_dp))
call assert(all(abs(lam(:1, 1) - [-2.280821_dp]) < 1e-6_dp))

end program
