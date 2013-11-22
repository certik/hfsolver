program test_sto_mg6

! Calculates the results for the STO article

use types, only: dp
use sto, only: stoints2, slater_sto_screen, sto_V_screen
use utils, only: assert
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
real(dp), parameter :: Dlist(*) = [1._dp, -1._dp]

Lmax = 2
allocate(nbfl(0:Lmax), nl(9, 0:Lmax), zl(9, 0:Lmax), focc(3, 0:Lmax))
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
ndof = sum(nbfl)
allocate(P_(n, n, 0:Lmax), C(n, n, 0:Lmax), H(n, n, 0:Lmax), lam(n, 0:Lmax))

open(newunit=u, file="mg.log", status="replace")
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
        call slater_sto_screen(nbfl, nl, zl, slater, D, &
            verbose=.true.)
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
