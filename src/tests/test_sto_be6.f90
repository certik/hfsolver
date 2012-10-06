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
integer :: n, Z, Nscf, Lmax, ndof, i, j, l, u
real(dp) :: alpha, Etot, tolE, tolP, Ekin
real(dp), allocatable :: H(:, :, :), P_(:, :, :), C(:, :, :), lam(:, :)
real(dp) :: D
!real(dp), parameter :: Dlist(*) = [1._dp, 10._dp, 100._dp, &
!    1e4_dp, 1e6_dp, -1._dp]
real(dp), parameter :: Dlist(*) = [1._dp, 1e6_dp, -1._dp]

Lmax = 2
allocate(nbfl(0:Lmax), nl(9, 0:Lmax), zl(9, 0:Lmax), focc(2, 0:Lmax))
nbfl = 0
focc = 0
focc(:2, 0) = [2, 2]
nbfl(0) = 5
nl(:5, 0) = [1, 1, 3, 2, 2]
zl(:5, 0) = [5.4297_dp, 2.9954_dp, 3.5810_dp, 1.1977_dp, 0.8923_dp]
nbfl(1) = 5
nl(:5, 1) = [2, 2, 4, 3, 3]
zl(:5, 1) = [5.6998_dp, 2.7850_dp, 4.1500_dp, 1.4387_dp, 0.9819_dp]
nbfl(2) = 2
nl(:2, 2) = [3, 3]
zl(:2, 2) = [1.2662_dp, 7.8314_dp]

Z = 4
tolE = 1e-10_dp
tolP = 1e-4_dp
alpha = 0.6_dp
Nscf = 100

n = maxval(nbfl)
ndof = sum(nbfl)
allocate(P_(n, n, 0:Lmax), C(n, n, 0:Lmax), H(n, n, 0:Lmax), lam(n, 0:Lmax))

open(newunit=u, file="be.log", status="replace")
do j = 1, size(Dlist)
    print *
    print *, "-----------------------------------------------------"
    print *, "total  DOFs =", ndof
    call stoints2(Z, nbfl, nl, zl, S, T, V, slater)
    D = Dlist(j)
    if (D < 0) then
        print *, "D = oo (Coulomb case)"
    else
        print *, "D =", D
        call sto_V_screen(Z, nbfl, nl, zl, V, D)
        call slater_sto_screen(nbfl, nl, zl, slater, D)
    end if

    H = T + V

    print *, "SCF cycle:"
    call doscf(nbfl, H, slater, S, focc, Nscf, tolE, tolP, alpha, C, P_, lam, Etot)
    Ekin = kinetic_energy(nbfl, P_, T)
!    call printall(nbfl, nl, zl, lam, C, Ekin, Etot)
    print *, "Debye length D =", D
    call printlam(nbfl, lam, Ekin, Etot)
    write(u, "(es23.16,' ',es23.16)") D, Etot
    do l = 0, ubound(nbfl, 1)
        do i = 1, nbfl(l)
            write(u, "(i4,i4,' ',es23.16)") i, l, lam(i, l)
        end do
    end do
end do
close(u)

end program
