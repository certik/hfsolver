module sto_mg_basis
use types, only: dp
use utils, only: assert
implicit none
private
public sto_optimized, sto_even_tempered, alpha_beta_step

contains

subroutine sto_optimized(Lmax, nbfl, nl, zl)
integer, intent(in) :: Lmax
integer, allocatable, intent(out) :: nl(:, :), nbfl(:)
real(dp), allocatable, intent(out) :: zl(:, :)
call assert(Lmax == 1)
allocate(nbfl(0:Lmax))
nbfl = [10, 7]
allocate(nl(maxval(nbfl), 0:Lmax), zl(maxval(nbfl), 0:Lmax))
nl(:nbfl(0), 0) = [2, 1, 2, 1, 1, 2, 1, 1, 2, 1]
zl(:nbfl(0), 0) = [35.988252_dp, 19.193253_dp, 16.093419_dp, 10.773377_dp, &
                7.853579_dp, 6.023976_dp, 2.905152_dp, 1.694174_dp, &
                0.750076_dp, 0.711410_dp]
nl(:nbfl(1), 1) = [3, 2, 3, 2, 2, 2, 2]
zl(:nbfl(1), 1) = [32.077221_dp, 12.495413_dp, 9.644163_dp, 5.042248_dp, &
                3.070882_dp, 2.086080_dp, 0.711410_dp]
end subroutine

subroutine sto_even_tempered(Lmax, nbfl, alpha, beta, nl, zl)
integer, intent(in) :: Lmax
integer, intent(in) :: nbfl(0:)
real(dp), intent(in) :: alpha(0:), beta(0:)
integer, allocatable, intent(out) :: nl(:, :)
real(dp), allocatable, intent(out) :: zl(:, :)
integer :: i, k
call assert(Lmax == 1)
allocate(nl(maxval(nbfl), 0:Lmax), zl(maxval(nbfl), 0:Lmax))
nl = 0
zl = 0
do i = 0, Lmax
    nl(:nbfl(i), i) = i+1
    zl(:nbfl(i), i) = [(alpha(i)*beta(i)**k, k=1,nbfl(i))]
end do
end subroutine

subroutine alpha_beta_step(alpha, beta, a, b, N)
real(dp), intent(inout) :: alpha, beta
real(dp), intent(in) :: a, b
integer, intent(in) :: N
real(dp) :: beta_old
call assert(a > 0)
call assert(-1 < b .and. b < 0)
beta_old = beta
beta = beta_old**((N/(N-1._dp))**b)
alpha = ((beta-1)/(beta_old-1))**a * alpha
end subroutine

end module


program sto_mg_conv

! Even-tempered (ET) STO basis convergence for Mg.

use types, only: dp
use sto, only: stoints2, get_basis2, slater_sto_screen, sto_V_screen
use utils, only: assert
use constants, only: pi, ang2bohr, Ha2eV
use radialscf, only: doscf, kinetic_energy, slater2int22, &
    get_basis, radialC2C, &
    radiallam2lam
use hfprint, only: printall, printlam
use mbpt, only: transform_int2, mbpt2, mbpt3, mbpt4, transform_int22
use gf2, only: find_poles, find_pole_diag, total_energy, plot_poles
use sorting, only: argsort
use scf, only: ijkl2intindex, ijkl2intindex2, create_intindex_sym4
use sto_mg_basis, only: sto_optimized, sto_even_tempered, alpha_beta_step
implicit none

integer, allocatable :: nl(:, :), nbfl(:)
real(dp), allocatable :: zl(:, :), focc(:, :)
real(dp), allocatable :: S(:, :, :), T(:, :, :), V(:, :, :), slater(:, :)
integer :: n, Z, m, Nscf, Lmax, ndof
real(dp) :: alpha_scf, Etot, tolE, tolP, Ekin
real(dp), allocatable :: H(:, :, :), P_(:, :, :), C(:, :, :), lam(:, :)
real(dp), allocatable :: alpha(:), beta(:)
integer :: Nb, u, i

Lmax = 1
allocate(focc(3, 0:Lmax))
focc = 0
focc(:3, 0) = [2, 2, 2]
focc(:1, 1) = [6]

Z = 12
tolE = 1e-10_dp
tolP = 1e-4_dp
alpha_scf = 0.6_dp
Nscf = 100

open(newunit=u, file="Etot.txt", status="replace")
close(u)

allocate(alpha(0:Lmax), beta(0:Lmax))
allocate(nbfl(0:Lmax))
alpha = [0.40938054_dp, 0.92139926_dp]
beta = [1.61259870_dp, 1.81117158_dp]
nbfl = [10, 5]

do Nb = 5, 30
    print *, "Nb =", Nb
    !call sto_optimized(Lmax, nbfl, nl, zl)
    call sto_even_tempered(Lmax, nbfl, alpha, beta, nl, zl)

    n = maxval(nbfl)
    ndof = sum(nbfl)
    print *, "total  DOFs =", ndof
    allocate(S(n, n, 0:Lmax), T(n, n, 0:Lmax), V(n, n, 0:Lmax))
    m = ndof*(ndof+1)/2
    allocate(slater(m*(m+1)/2, 0:2*Lmax))
    call stoints2(Z, nbfl, nl, zl, S, T, V, slater)
    allocate(P_(n, n, 0:Lmax), C(n, n, 0:Lmax), H(n, n, 0:Lmax), lam(n, 0:Lmax))

    H = T + V
    print *, "SCF cycle:"
    call doscf(nbfl, H, slater, S, focc, Nscf, tolE, tolP, alpha_scf, C, P_, lam, Etot)
    Ekin = kinetic_energy(nbfl, P_, T)
    !call printall(nbfl, nl, zl, lam, C, Ekin, Etot)
    call printlam(nbfl, lam, Ekin, Etot)

    open(newunit=u, file="Etot.txt", position="append", status="old")
    write(u, "(i4, ' ', i5, ' ', es23.16, ' ', es23.16, ' ', 2f10.4, 2f10.4)") Nb, ndof, Etot, &
        Ekin+Etot, alpha, beta
    close(u)

    do i = 1, 2
        nbfl(0) = nbfl(0) + 1
        call alpha_beta_step(alpha(0), beta(0), 0.5_dp, -0.5_dp, nbfl(0))
    end do
    do i = 1, 1
        nbfl(1) = nbfl(1) + 1
        call alpha_beta_step(alpha(1), beta(1), 0.6_dp, -0.45_dp, nbfl(1))
    end do

    deallocate(S, T, V, slater, P_, C, H, lam)
end do

end program
