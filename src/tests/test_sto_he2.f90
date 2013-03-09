program test_sto_he2
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

allocate(nbfl(0:0), nl(9, 0:0), zl(9, 0:0), focc(1, 0:0))
nbfl = 0
focc = 0
focc(:1, 0) = [2]
nbfl(0) = 9
nl(:9, 0) = [1, 1, 1, 2, 2, 2, 3, 3, 3]
zl(:9, 0) = [1.4_dp, 2.5_dp, 0.67_dp, 4.3_dp, 1._dp, 0.67_dp, 5.5_dp, &
    3.9_dp, 0.67_dp]

Z = 2
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

call assert(abs(Etot - (-2.86167998_dp)) < 1e-8_dp)
call assert(abs(Ekin - (+2.86167931_dp)) < 1e-8_dp)
call assert(abs(lam(1, 0) - (-0.917956_dp)) < 1e-6_dp)

end program
