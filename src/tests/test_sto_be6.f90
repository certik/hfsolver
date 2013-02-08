program test_sto_be6

! Calculates the results for the STO article

use types, only: dp
use sto, only: stoints2, slater_sto_screen, sto_V_screen
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

real(dp), allocatable :: S(:, :, :), T(:, :, :), V(:, :, :), slater(:, :)
integer :: n, Z, Nscf, Lmax, ndof, j, u
real(dp) :: alpha, Etot, tolE, tolP, Ekin
real(dp), allocatable :: H(:, :, :), P_(:, :, :), C(:, :, :), lam(:, :)
real(dp) :: D
real(dp), parameter :: Dlist(*) = [-1._dp, 20._dp, 10._dp, 8._dp, 6._dp, &
    5._dp, 4._dp, 3._dp, 2._dp, 1.5_dp, 1.25_dp, 1.15_dp, 1.05_dp]

Lmax = 0
allocate(nbfl(0:Lmax), nl(9, 0:Lmax), zl(9, 0:Lmax), focc(2, 0:Lmax))
nbfl = 0
focc = 0
focc(:2, 0) = [2]
nbfl(0) = 8
nl(:8, 0) = [2, 1, 2, 1, 1, 2, 1, 1]
zl(:8, 0) = [18.890445_dp, 9.238787_dp, 7.517513_dp, 5.100368_dp, &
    3.276630_dp, 2.270243_dp, 1.192963_dp, 0.930957_dp]

Z = 6
tolE = 1e-10_dp
tolP = 1e-4_dp
alpha = 0.6_dp
Nscf = 100

n = maxval(nbfl)
ndof = sum(nbfl)
allocate(P_(n, n, 0:Lmax), C(n, n, 0:Lmax), H(n, n, 0:Lmax), lam(n, 0:Lmax))

open(newunit=u, file="be.log", status="replace")
write(u, "(a10,' ',a10,' ',a12)") "D", "mu", "Model A"
do j = 1, size(Dlist)
    print *
    print *, "-----------------------------------------------------"
    print *, "total  DOFs =", ndof
    call stoints2(Z, nbfl, nl, zl, S, T, V, slater)
    D = Dlist(j)
    if (D < 0) then
        print *, "D = oo (Coulomb case)"
    else
        print *, "D =", D, "lam = 1/D =", 1/D
        call sto_V_screen(Z, nbfl, nl, zl, V, D)
!        call slater_sto_screen(nbfl, nl, zl, slater, D)
    end if

    H = T + V

    print *, "SCF cycle:"
    call doscf(nbfl, H, slater, S, focc, Nscf, tolE, tolP, alpha, C, P_, lam, Etot)
    Ekin = kinetic_energy(nbfl, P_, T)
!    call printall(nbfl, nl, zl, lam, C, Ekin, Etot)
    print *, "Debye length D =", D
    call printlam(nbfl, lam, Ekin, Etot)
    write(u, "(f10.2,' ',f10.3,' ',f12.4)") D, 1/D, Etot
end do
close(u)

end program
