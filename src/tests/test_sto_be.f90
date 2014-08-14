program test_sto_be
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

allocate(nbfl(0:2), nl(12, 0:2), zl(12, 0:2), focc(2, 0:2))
!allocate(nbfl(0:3), nl(12, 0:3), zl(12, 0:3), focc(2, 0:3))
nbfl = 0
focc = 0
focc(:2, 0) = [2, 2]
nbfl(0) = 12
nl(:12, 0) = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3]
zl(:12, 0) = [0.1_dp, 1.5_dp, 3.4703_dp, 6.3681_dp, &
    0.7516_dp, 0.9084_dp, 5.80765_dp, 8.53828_dp, &
    0.08_dp, 16.1268_dp, 20.2055_dp, 31.4944_dp]
nbfl(1) = 12
nl(:12, 1) = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]
zl(:12, 1) = [0.20135_dp, 0.95703_dp, 1.943_dp, 2.6986_dp, &
    0.09_dp, 1.46565_dp, 3.35114_dp, 5.1389_dp, &
    0.15969_dp, 0.75902_dp, 8.80274_dp, 11.521_dp]
nbfl(2) = 12
nl(:12, 2) = [3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5]
zl(:12, 2) = [0.074284_dp, 0.373267_dp, 0.83752_dp, 1.3652_dp, &
    1.8273_dp, 2.30093_dp, 0.0805239_dp, 0.512334_dp, &
    3.22431_dp, 4.47682_dp, 6.2468_dp, 8.68472_dp]
!nbfl(3) = 12
!nl(:12, 3) = [4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5]
!zl(:12, 3) = [14.64884_dp, 1.17503_dp, 17.49577_dp, 4.91401_dp, 9.22804_dp, &
!    11.94673_dp, 2.73331_dp, 3.15699_dp, 5.23396_dp, 0.71385_dp, 20.02434_dp, &
!    0.84142_dp]
!zl = zl * 0.4_dp

Z = 4
tolE = 1e-10_dp
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
call printall(nbfl, nl, zl, lam, C, Ekin, Etot)
call printlam(nbfl, lam, Ekin, Etot)

call assert(abs(Etot - (-14.57279688_dp)) < 1e-8_dp)
call assert(all(abs(lam(:2, 0) - [-4.733142_dp, -0.309212_dp]) < 1e-6_dp))
call assert(all(abs(lam(10:, 1) - [18.840697_dp, 47.564325_dp, &
    146.358091_dp]) < 1e-6_dp))
call assert(all(abs(lam(10:, 2) - [8.393357_dp, 16.676739_dp, &
    31.940327_dp]) < 1e-6_dp))

end program
