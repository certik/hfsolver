program sto_be_screen_gf

! Screened Green's function calculation
! The results for D->oo are to be compared against [1].

! [1] Doll, J. D., & Reinhardt, W. P. (1972). Many-Body Green’s Functions for
! Finite, Nonuniform Systems: Applications to Closed Shell Atoms. The Journal
! of Chemical Physics, 57(3), 1169–1184. doi:10.1063/1.1678374

use types, only: dp
use sto, only: stoints2, get_basis2, slater_sto_screen, sto_V_screen
use utils, only: assert
use radialscf, only: doscf, kinetic_energy, slater2int22, &
    get_basis, radialC2C, &
    radiallam2lam
use hfprint, only: printall, printlam
use mbpt, only: transform_int2, mbpt2, mbpt3, mbpt4, transform_int22
use gf2, only: find_poles, find_pole_diag, total_energy, plot_poles
use sorting, only: argsort
use scf, only: ijkl2intindex, ijkl2intindex2, create_intindex_sym4
implicit none

integer, allocatable :: nl(:, :), nbfl(:)
real(dp), allocatable :: zl(:, :), focc(:, :)

real(dp), allocatable :: S(:, :, :), T(:, :, :), V(:, :, :), slater(:, :)
integer :: n, Z, m, Nscf, Lmax, ndof
real(dp) :: alpha, Etot, tolE, tolP, Ekin
real(dp), allocatable :: H(:, :, :), P_(:, :, :), C(:, :, :), lam(:, :)

integer, allocatable :: nlist(:), llist(:), mlist(:)
real(dp), allocatable :: int2(:), moint2(:), Ctot(:, :), lamtot(:), Egreen(:)
integer, allocatable :: intindex(:, :, :, :)
integer :: Nelec

real(dp) :: D
integer, allocatable :: idx(:)
integer :: i, it, it_N, u, N_elec
real(dp) :: Dmin, Dmax, it_a, it_b

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
N_elec = round(sum(focc))
print *, "total  DOFs =", ndof
allocate(S(n, n, 0:Lmax), T(n, n, 0:Lmax), V(n, n, 0:Lmax))
m = ndof*(ndof+1)/2
allocate(slater(m*(m+1)/2, 0:2*Lmax))
call stoints2(Z, nbfl, nl, zl, S, T, V, slater)
allocate(P_(n, n, 0:Lmax), C(n, n, 0:Lmax), H(n, n, 0:Lmax), lam(n, 0:Lmax))
call get_basis(nbfl, nlist, llist, mlist)
allocate(Ctot(size(nlist), size(nlist)), lamtot(size(nlist)), idx(size(nlist)))
allocate(Egreen(N_elec))
m = size(nlist)**2 * (size(nlist)**2 + 3) / 4
allocate(int2(m), moint2(m))
allocate(intindex(size(nlist), size(nlist), size(nlist), size(nlist)))
call create_intindex_sym4(intindex)

open(newunit=u, file="eigs_D.txt", status="replace")
close(u)
Dmin = 1._dp
Dmax = 1e5
it_N = 100
it_a = Dmin * (Dmin/Dmax)**(1._dp/(it_N-1))
it_b = log(Dmax/Dmin)/(it_N-1)
do it = 1, it_N
    D = it_a*exp(it_b*it)
    print *, "D =", D
    call sto_V_screen(Z, nbfl, nl, zl, V, D)
    call slater_sto_screen(nbfl, nl, zl, slater, D)
    H = T + V
    print *, "SCF cycle:"
    call doscf(nbfl, H, slater, S, focc, Nscf, tolE, tolP, alpha, C, P_, lam, Etot)
    Ekin = kinetic_energy(nbfl, P_, T)
!    call printall(nbfl, nl, zl, lam, C, Ekin, Etot)
    call printlam(nbfl, lam, Ekin, Etot)
    call radiallam2lam(nlist, llist, lam, lamtot)
    idx = argsort(lamtot)
    nlist = nlist(idx)
    llist = llist(idx)
    mlist = mlist(idx)
    lamtot = lamtot(idx)
    call slater2int22(nbfl, nlist, llist, mlist, slater, int2)
    call radialC2C(nlist, llist, mlist, C, Ctot)
    moint2 = transform_int22(int2, intindex, Ctot)
    Nelec = Z
    print *, "Green's function calculation:"
    print *, "i         lam(i)            dE         lam(i)+dE"
    Egreen = 0
    do i = 1, size(Egreen)
        Egreen(i) = find_pole_diag(i, moint2, intindex, lamtot, Nelec/2, &
            200, 1e-10_dp)
        print "(i4, f15.6, f15.6, f15.6)", i, lamtot(i), Egreen(i)-lamtot(i), &
            Egreen(i)
    end do

    open(newunit=u, file="be_eigs_D.txt", position="append", status="old")
    write(u, "(6(es23.16, ' '))") D, Etot, lamtot(:2), Egreen(:2)
    close(u)

end do

contains

integer function round(x) result(r)
! Rounds "x" to the nearest integer
real(dp), intent(in) :: x
r = int(x + 0.5_dp)
end function

end program
