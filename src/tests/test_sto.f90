program test_sto
use types, only: dp
use sto, only: get_basis, stoints
use utils, only: assert
use constants, only: Ha2eV
use scf, only: doscf, get_nuclear_energy
implicit none

integer, allocatable :: nl(:), ll(:), ml(:)
real(dp), allocatable :: zl(:)

real(dp), allocatable :: S(:, :), T(:, :), V(:, :), int2(:)
integer :: n, Z, m, i, Nelec, Nscf
real(dp) :: alpha, Etot, tolE, tolP, Ekin
real(dp), allocatable :: H(:, :), P_(:, :), C(:, :), lam(:)

allocate(nl(40), zl(40), ll(40), ml(40))
ml = 0
ll(8*0+1:8*1) = 0 ! s
nl(8*0+1:8*1) = [1, 2, 2, 3, 3, 3, 3, 4]
zl(8*0+1:8*1) = [12.0_dp, 13.6_dp, 9.3_dp, 6.50_dp, 4.2_dp, 1.4_dp, 0.9_dp, 2.5_dp]
ll(8*1+1:8*2) = 1 ! p 3x
nl(8*1+1:8*2) = [2, 2, 2, 3, 3, 3, 3, 4]
zl(8*1+1:8*2) = [12.5_dp,  9.2_dp, 6.0_dp, 0.50_dp, 5.3_dp, 3.7_dp, 2.5_dp, 1.2_dp]
ml(8*1+1:8*2) = -1
ll(8*2+1:8*3) = ll(8*1+1:8*2)
nl(8*2+1:8*3) = nl(8*1+1:8*2)
zl(8*2+1:8*3) = zl(8*1+1:8*2)
ml(8*2+1:8*3) =  0
ll(8*3+1:8*4) = ll(8*1+1:8*2)
nl(8*3+1:8*4) = nl(8*1+1:8*2)
zl(8*3+1:8*4) = zl(8*1+1:8*2)
ml(8*3+1:8*4) = +1
ll(8*4+1:8*5) = 2 ! d
nl(8*4+1:8*5) = [3, 3, 3, 3, 4, 4, 4, 4]
zl(8*4+1:8*5) = [12.8_dp,  9.6_dp, 6.4_dp, 0.75_dp, 5.2_dp, 3.6_dp, 2.1_dp, 1.4_dp]

Z = 12
Nelec = 12
tolE = 1e-6_dp
tolP = 1e-4_dp
alpha = 0.6_dp
Nscf = 100

n = size(nl)
print *, "total  DOFs =", n
allocate(S(n, n), T(n, n), V(n, n))
m = n*(n+1)/2
allocate(int2(m*(m+1)/2))
call stoints(Z, nl, zl, ll, ml, S, T, V, int2)

allocate(P_(n, n), C(n, n), lam(n))


H = T + V
!print "(8f5.2)", T(:8, :8)
!print *
!print "(8f5.2)", V(:8, :8)
!stop "O"


!call print_STH()
!call print_int2()
!stop

print *, "SCF cycle:"
call doscf(H, int2, S, Nelec/2, Nscf, tolE, tolP, alpha, C, P_, lam, Etot)
Ekin = sum(P_*T)
print *, "Orbital energies:"
print *, "  i          E [a.u.]          E [eV]"
do i = 1, size(lam)
    print "(i4,f18.6,f18.6)", i, lam(i), lam(i) * Ha2eV
end do
print "(a,e16.8)", "        EKIN+EHF (a.u.):", Ekin + Etot
print "(a,e16.8)", "  KINETIC ENERGY (a.u.):", Ekin
print "(a,e16.8)", "HF ATOMIC ENERGY (a.u.):", Etot

contains

subroutine print_STH()
integer :: i, j
do j = 1, size(S, 1)
    do i = j, size(S, 1)
        if (ll(i) == ll(j) .and. ml(i) == ml(j) .and. ml(i) == 0) then
            print "(3e25.13)", S(i, j), T(i, j), H(i, j)
        end if
    end do
end do
end subroutine

subroutine print_int2()
integer :: i, u
open(newunit=u, file="int2", status="replace")
do i = 1, size(int2)
    if (abs(int2(i)) > tiny(1._dp)) write(u,"(e25.13)") int2(i)
end do
close(u)
end subroutine

end program
