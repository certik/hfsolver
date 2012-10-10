program sto_be_screen_gf

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
use scf, only: ijkl2intindex, ijkl2intindex2
implicit none

integer, allocatable :: nl(:, :), nbfl(:)
real(dp), allocatable :: zl(:, :), focc(:, :)

real(dp), allocatable :: S(:, :, :), T(:, :, :), V(:, :, :), slater(:, :)
integer :: n, Z, m, Nscf, Lmax, ndof
real(dp) :: alpha, Etot, tolE, tolP, Ekin
real(dp), allocatable :: H(:, :, :), P_(:, :, :), C(:, :, :), lam(:, :)

integer, allocatable :: nlist(:), llist(:), mlist(:)
real(dp), allocatable :: int2(:), moint2(:), Ctot(:, :), Htot(:, :), lamtot(:)
integer, allocatable :: intindex(:, :, :, :)
integer :: Nelec

!real(dp) :: E2, E3, E4
real(dp) :: Ntot, Egreen, D
!real(dp), allocatable :: lam_green(:)
integer, allocatable :: idx(:)
integer :: i

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
print *, "total  DOFs =", ndof
allocate(S(n, n, 0:Lmax), T(n, n, 0:Lmax), V(n, n, 0:Lmax))
m = ndof*(ndof+1)/2
allocate(slater(m*(m+1)/2, 0:2*Lmax))
call stoints2(Z, nbfl, nl, zl, S, T, V, slater)

D = 20
print *, "D =", D
call sto_V_screen(Z, nbfl, nl, zl, V, D)
call slater_sto_screen(nbfl, nl, zl, slater, D)

allocate(P_(n, n, 0:Lmax), C(n, n, 0:Lmax), H(n, n, 0:Lmax), lam(n, 0:Lmax))

H = T + V

print *, "SCF cycle:"
call doscf(nbfl, H, slater, S, focc, Nscf, tolE, tolP, alpha, C, P_, lam, Etot)
Ekin = kinetic_energy(nbfl, P_, T)
call printall(nbfl, nl, zl, lam, C, Ekin, Etot)
call printlam(nbfl, lam, Ekin, Etot)

!call assert(abs(Etot - (-14.57278856_dp)) < 1e-8_dp)

call get_basis(nbfl, nlist, llist, mlist)
allocate(Ctot(size(nlist), size(nlist)), lamtot(size(nlist)), idx(size(nlist)))
allocate(Htot(size(nlist), size(nlist)))
call radiallam2lam(nlist, llist, lam, lamtot)

m = size(nlist)**2 * (size(nlist)**2 + 3) / 4
allocate(int2(m), moint2(m))
allocate(intindex(size(nlist), size(nlist), size(nlist), size(nlist)))
call create_index(intindex)


! permute() is not needed any more:
!call permute(int2, idx, size(nlist))
idx = argsort(lamtot)
nlist = nlist(idx)
llist = llist(idx)
mlist = mlist(idx)
lamtot = lamtot(idx)
call slater2int22(nbfl, nlist, llist, mlist, slater, int2)

call radialC2C(nlist, llist, mlist, C, Ctot)
call radialC2C(nlist, llist, mlist, H, Htot)

Htot = transform_H(Htot, Ctot)

moint2 = transform_int22(int2, intindex, Ctot)

Nelec = Z
!print *, "Calculating MPBT 2"
!E2 = mbpt2(moint2, lamtot, Nelec/2)
!print *, "Calculating MPBT 3"
!E3 = mbpt3(moint2, lamtot, Nelec/2)
!print *, "Calculating MPBT 4"
!E4 = mbpt4(moint2, lamtot, Nelec/2)

print *, "MBPT results:"
print "(' E0+E1 (HF)    = ',f15.8)", Etot
!print "(' E2    (MBPT2) = ',f15.8)", E2
!print "(' E3    (MBPT3) = ',f15.8)", E3
!print "(' E4    (MBPT4) = ',f15.8)", E4
!print "(' E0+E1+E2      = ',f15.8)", Etot + E2
!print "(' E0+E1+E2+E3   = ',f15.8)", Etot + E2 + E3
!print "(' E0+E1+E2+E3+E4= ',f15.8)", Etot + E2 + E3 + E4
!call assert(abs(E2 - (-0.03532286_dp)) < 1e-8_dp)

print *, "Green's function calculation:"
print *, "i         lam(i)            dE         lam(i)+dE"
do i = 1, size(lam)
    if (lamtot(i) > 0) exit
    Egreen = find_pole_diag(i, moint2, intindex, lamtot, Nelec/2, 200, 1e-10_dp)
    print "(i4, f15.6, f15.6, f15.6)", i, lamtot(i), Egreen-lamtot(i), Egreen
    if (i == 1) then
        ! [1] has energy:         -4.612
        !call assert(abs(Egreen - (-4.61175302_dp)) < 1e-8_dp)
    else if (i == 2) then
        ! [1] has energy:         -0.327
        !call assert(abs(Egreen - (-0.32724971_dp)) < 1e-8_dp)
    end if
end do

!call plot_poles(moint2, intindex, lamtot, Nelec/2)

print *, "GF total energy (p=20):"
call total_energy(moint2, intindex, lamtot, Nelec/2, 20, Htot, &
        200, 0.1_dp, 0.1_dp, 0._dp, Ntot, Egreen)
print *, "Ntot =", Ntot
print *, "Etot =", Egreen
!do i = 5, 60, 5
!    print *, "p =", i
!    call total_energy(moint2, intindex, lamtot, Nelec/2, i, Htot, &
!        200, 0.1_dp, 0.1_dp, 0._dp, Ntot, Egreen)
!    print *, "Ntot =", Ntot
!    print *, "Etot =", Egreen
!end do
! These numbers are converged to at least 1e-9:
!call assert(abs(Ntot   - (  4.004281307_dp)) < 1e-6_dp)
!call assert(abs(Egreen - (-14.650438751_dp)) < 1e-6_dp)
! [1] has total energy:   -14.641639
! Note: if we extended the countour to -55 a.u., we get:
! In [1]: -14.666238504287749 / 4.0052366660634684 * 4
! Out[1]: -14.6470630597741
! So after normalization by the number of particles, we get a close energy.


!allocate(lam_green(size(nlist)))
!lam_green = find_poles(moint2, lamtot, Nelec/2, 20)
!do i = 1, size(lamtot)
!    if (lamtot(i) > 0) exit
!    print "(i4, f15.6, f15.6, f15.6)", i, lamtot(i), lam_green(i)-lamtot(i), &
!        lam_green(i)
!end do

contains


subroutine permute(A, idx, n)
real(dp), intent(inout) :: A(:)
integer, intent(in) :: idx(:), n
real(dp) :: B(size(A))
integer :: i, j, k, l
do i = 1, n
    do j = 1, i
        do k = 1, n
            do l = 1, k
                if ((i-1)*i/2 + j < (k-1)*k/2 + l) cycle
                B(ijkl2intindex(i, j, k, l)) = &
                    A(ijkl2intindex(idx(i), idx(j), idx(k), idx(l)))
            end do
        end do
    end do
end do
A = B
end subroutine

function transform_H(H, C) result(Horb)
real(dp), intent(in) :: H(:, :), C(:, :)
real(dp) :: Horb(size(H, 1), size(H, 2))
Horb = matmul(matmul(transpose(C), H), C)
end function

subroutine create_index(intindex)
integer, intent(out) :: intindex(:, :, :, :)
integer :: i, j, k, l, ijkl, n
n = size(intindex, 1)
! Uncomment the commented lines to run checks:
!intindex = -1
ijkl = 1
do i = 1, n
    do j = 1, i
        do k = 1, i
            do l = 1, i
                if (i == j .and. k < l) cycle
                if (i == k .and. j < l) cycle
                if (i == l .and. j < k) cycle
!                call assert(ijkl2intindex2(i, j, k, l) == ijkl)
                intindex(i, j, k, l) = ijkl
                intindex(j, i, l, k) = ijkl
                intindex(k, l, i, j) = ijkl
                intindex(l, k, j, i) = ijkl
                ijkl = ijkl + 1
            end do
        end do
    end do
end do
!call assert(.not. (any(intindex == -1)))
end subroutine

end program
