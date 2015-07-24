program test_sto_hm
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

allocate(nbfl(0:2), nl(5, 0:2), zl(5, 0:2), focc(2, 0:2))
nbfl = 0
focc = 0
focc(1, 0) = 2
nbfl(0) = 5
nl(:5, 0) = [1, 1, 2, 3, 3]
zl(:5, 0) = [0.28382_dp, 0.51444_dp, 0.8525_dp, 0.79958_dp, 1.09728_dp]
nbfl(1) = 4
nl(:4, 1) = [2, 2, 3, 4]
zl(:4, 1) = [0.103336_dp, 0.145652_dp, 0.221232_dp, 0.228868_dp] * 5
nbfl(2) = 3
nl(:3, 2) = [3, 3, 4]
zl(:3, 2) = [0.14546_dp, 0.193412_dp, 0.278776_dp] * 5

Z = 1
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
call printall(nbfl, nl, zl, lam, C, Ekin, Etot)
call printlam(nbfl, lam, Ekin, Etot)

call assert(abs(Etot - (-0.48789731_dp)) < 1e-8_dp)
call assert(all(abs(lam(:, 0) - [-0.046227_dp, 0.126799_dp, 0.342570_dp, &
    0.914781_dp, 4.111577_dp]) < 1e-6_dp))
call assert(all(abs(lam(:4, 1) - [0.186112_dp, 0.422181_dp, 0.964922_dp, &
    2.399488_dp]) < 1e-6_dp))
call assert(all(abs(lam(:3, 2) - [0.325247_dp, 0.966288_dp, &
    2.656610_dp]) < 1e-6_dp))

end program
