program test_sto_ra

! The STO basis is taken from (and the energies compared against):
!
! THEORETICAL CHEMISTRY ACCOUNTS: THEORY, COMPUTATION, AND MODELING (THEORETICA
! CHIMICA ACTA)
! Volume 104, Number 5 (2000), 411-413, DOI: 10.1007/s002140000150
!
! http://www.springerlink.com/content/v7p00qyqgh5t3neu/
! http://www.unb.ca/fredericton/science/chem/ajit/download.htm

use types, only: dp
use sto, only: stoints2, get_basis2, get_values
use utils, only: assert, savetxt
use mesh, only: linspace
use radialscf, only: doscf, kinetic_energy
use hfprint, only: printall, printlam
implicit none

integer, allocatable :: nl(:, :), nbfl(:)
real(dp), allocatable :: zl(:, :), focc(:, :)

real(dp), allocatable :: S(:, :, :), T(:, :, :), V(:, :, :), slater(:, :)
integer :: n, Z, m, Nscf, Lmax, ndof
real(dp) :: alpha, Etot, tolE, tolP, Ekin
real(dp), allocatable :: H(:, :, :), P_(:, :, :), C(:, :, :), lam(:, :)
real(dp), allocatable :: R(:), u(:, :, :)

allocate(nbfl(0:3), nl(15, 0:3), zl(15, 0:3), focc(7, 0:3))
nbfl = 0
focc = 0
focc(:7, 0) = [2, 2, 2, 2, 2, 2, 2]
focc(:5, 1) = [6, 6, 6, 6, 6]
focc(:3, 2) = [10, 10, 10]
focc(:1, 3) = [14]
nbfl(0) = 15
nl(:15, 0) = [1, 1, 2, 2, 3, 3, 3, 4, 3, 3, 3, 3, 3, 3, 3]
zl(:15, 0) = [102.772371_dp, 85.593215_dp, 70.820905_dp, 41.242337_dp, &
    36.905364_dp, 25.990833_dp, 20.824259_dp, 10.834095_dp, 10.835315_dp, &
    5.559242_dp, 4.498561_dp, 2.306842_dp, 1.710801_dp, 0.761325_dp, &
    0.569949_dp]
nbfl(1) = 12
nl(:12, 1) = [2, 2, 3, 3, 4, 4, 3, 3, 3, 3, 3, 3]
zl(:12, 1) = [83.018450_dp, 46.836527_dp, 41.801671_dp, 25.187270_dp, &
    25.102935_dp, 13.925446_dp, 13.918592_dp, 9.852338_dp, 5.415789_dp, &
    4.328539_dp, 1.999556_dp, 1.359217_dp]
nbfl(2) = 9
nl(:9, 2) = [3, 3, 4, 4, 4, 4, 4, 4, 4]
zl(:9, 2) = [54.653258_dp, 31.812073_dp, 26.239889_dp, 17.292475_dp, &
    13.510169_dp, 10.453460_dp, 7.146561_dp, 4.838823_dp, 3.220235_dp]
nbfl(3) = 6
nl(:6, 3) = [4, 4, 4, 5, 4, 4]
zl(:6, 3) = [30.150869_dp, 18.531819_dp, 12.237536_dp, 8.036556_dp, &
    7.428155_dp, 4.881720_dp]

Z = 88
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
allocate(R(20000))
allocate(u(size(R), size(C, 2), 0:ubound(C, 3)))
! This calculates the values of eigenvectors (for plotting)
!R = linspace(0._dp, 30._dp, size(R))
!u = get_values(nbfl, nl, zl, C, R)
!open(newunit=n, file="R.txt", status="replace")
!write(n, *) R
!close(n)
!call savetxt("eigs0.txt", u(:, :nbfl(0), 0))
!call savetxt("eigs1.txt", u(:, :nbfl(1), 1))
!call savetxt("eigs2.txt", u(:, :nbfl(2), 2))
!call savetxt("eigs3.txt", u(:, :nbfl(2), 3))
Ekin = kinetic_energy(nbfl, P_, T)
call printall(nbfl, nl, zl, lam, C, Ekin, Etot)
call printlam(nbfl, lam, Ekin, Etot)

call assert(abs(Etot - (-23094.30349262_dp)) < 1e-8_dp)
call assert(abs(Ekin - (+23094.30362539_dp)) < 1e-8_dp)
call assert(all(abs(lam(:7, 0) - [-3388.941040_dp, -587.744207_dp, &
    -147.871730_dp, -37.341571_dp, -8.253170_dp, -1.370676_dp, &
    -0.148760_dp]) < 1e-6_dp))
call assert(all(abs(lam(:1, 3) - [-12.430502_dp]) < 1e-6_dp))

end program
