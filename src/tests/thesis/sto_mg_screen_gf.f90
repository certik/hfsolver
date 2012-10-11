program sto_mg_screen_gf

! Screened Green's function calculation
! The results for D->oo are to be compared against [1].

! [1] Doll, J. D., & Reinhardt, W. P. (1972). Many-Body Green’s Functions for
! Finite, Nonuniform Systems: Applications to Closed Shell Atoms. The Journal
! of Chemical Physics, 57(3), 1169–1184. doi:10.1063/1.1678374

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
integer :: i, it, it_N, u
real(dp) :: Dmin, Dmax, it_a, it_b, t1, t2

Lmax = 1
allocate(nbfl(0:Lmax))
nbfl = [10, 7]
allocate(nl(maxval(nbfl), 0:Lmax), zl(maxval(nbfl), 0:Lmax), focc(3, 0:Lmax))
focc = 0
focc(:3, 0) = [2, 2, 2]
focc(:1, 1) = [6]
nl(:nbfl(0), 0) = [2, 1, 2, 1, 1, 2, 1, 1, 2, 1]
zl(:nbfl(0), 0) = [35.988252_dp, 19.193253_dp, 16.093419_dp, 10.773377_dp, &
                7.853579_dp, 6.023976_dp, 2.905152_dp, 1.694174_dp, &
                0.750076_dp, 0.711410_dp]
nl(:nbfl(1), 1) = [3, 2, 3, 2, 2, 2, 2]
zl(:nbfl(1), 1) = [32.077221_dp, 12.495413_dp, 9.644163_dp, 5.042248_dp, &
                3.070882_dp, 2.086080_dp, 0.711410_dp]

Z = 12
tolE = 1e-10_dp
tolP = 1e-4_dp
alpha = 0.6_dp
Nscf = 100

n = maxval(nbfl)
ndof = sum(nbfl)
Nelec = round(sum(focc))
print *, "total  DOFs =", ndof
allocate(S(n, n, 0:Lmax), T(n, n, 0:Lmax), V(n, n, 0:Lmax))
m = ndof*(ndof+1)/2
allocate(slater(m*(m+1)/2, 0:2*Lmax))
call stoints2(Z, nbfl, nl, zl, S, T, V, slater)
allocate(P_(n, n, 0:Lmax), C(n, n, 0:Lmax), H(n, n, 0:Lmax), lam(n, 0:Lmax))
call get_basis(nbfl, nlist, llist, mlist)
allocate(Ctot(size(nlist), size(nlist)), lamtot(size(nlist)), idx(size(nlist)))
allocate(Egreen(Nelec))
m = size(nlist)**2 * (size(nlist)**2 + 3) / 4
allocate(int2(m), moint2(m))
allocate(intindex(size(nlist), size(nlist), size(nlist), size(nlist)))
call create_intindex_sym4(intindex)

open(newunit=u, file="mg_eigs_D.txt", status="replace")
close(u)
Dmin = 1._dp
Dmax = 1e5
it_N = 100
it_a = Dmin * (Dmin/Dmax)**(1._dp/(it_N-1))
it_b = log(Dmax/Dmin)/(it_N-1)
do it = it_N, 1, -1
    D = it_a*exp(it_b*it)
    print *, "D =", D
    call sto_V_screen(Z, nbfl, nl, zl, V, D)
    call slater_sto_screen(nbfl, nl, zl, slater, D)
    H = T + V
    print *, "SCF cycle:"
    call doscf(nbfl, H, slater, S, focc, Nscf, tolE, tolP, alpha, C, P_, lam, Etot)
    Ekin = kinetic_energy(nbfl, P_, T)
    !call printall(nbfl, nl, zl, lam, C, Ekin, Etot)
    call printlam(nbfl, lam, Ekin, Etot)
    call radiallam2lam(nlist, llist, lam, lamtot)
    idx = argsort(lamtot)
    nlist = nlist(idx)
    llist = llist(idx)
    mlist = mlist(idx)
    lamtot = lamtot(idx)
    call slater2int22(nbfl, nlist, llist, mlist, slater, int2)
    call radialC2C(nlist, llist, mlist, C, Ctot)
    call cpu_time(t1)
    moint2 = transform_int22(int2, intindex, Ctot)
    call cpu_time(t2)
    print *, t2-t1
!    stop "OK"
    print *, "Green's function calculation:"
    print *, "i         lam(i)            dE         lam(i)+dE"
    Egreen = 0
    do i = 1, Nelec
        Egreen(i) = find_pole_diag(i, moint2, intindex, lamtot, Nelec/2, &
            200, 1e-10_dp)
        print "(i4, f15.6, f15.6, f15.6)", i, lamtot(i), Egreen(i)-lamtot(i), &
            Egreen(i)
    end do

    open(newunit=u, file="mg_eigs_D.txt", position="append", status="old")
    write(u, "(100(es23.16, ' '))") D, Etot, lamtot(:Nelec), Egreen(:Nelec)
    close(u)

end do

contains

integer function round(x) result(r)
! Rounds "x" to the nearest integer
real(dp), intent(in) :: x
r = int(x + 0.5_dp)
end function

end program
