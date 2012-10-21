module linalg_feast
use types, only: dp
use utils, only: stop_error, assert
use feast, only: feastinit, dfeast_sygv, dfeast_syev
implicit none
private
public eigh

interface eigh
    module procedure deigh_generalized
    module procedure deigh_simple
end interface eigh


contains

subroutine deigh_generalized(A, B, Emin, Emax, M0, lam, c)
real(dp), intent(in) :: A(:, :), B(:, :), Emin, Emax
integer, intent(in) :: M0
real(dp), intent(out) :: lam(:), c(:, :)
integer :: feastparam(64)
real(dp), allocatable :: res(:)
real(dp) :: epsout
integer :: i, loop, info, M, N, LDA
N = size(A, 1)
LDA = N
allocate(res(M0))

call feastinit(feastparam)
feastparam(1)=1 !! change from default value
call dfeast_sygv('F',N,A,LDA,B,LDA,feastparam,epsout,loop,Emin,Emax,M0,lam,c, &
    M,res,info)

if (info /= 0) then
    print *,'FEAST OUTPUT INFO', info
    call stop_error("info /= 0")
end if

print *, '# Search interval [Emin,Emax]', Emin, Emax
print *, '# mode found/subspace', M, M0
print *, '# iterations', loop
print *, 'TRACE', sum(lam(1:M))
print *, 'Relative error on the Trace', epsout
print *, 'Eigenvalues/Residuals'
do i = 1, M
    print *, i, lam(i), res(i)
end do

if (M < M0) call stop_error("M < M0")

end subroutine


subroutine deigh_simple(A, Emin, Emax, M0, lam, c)
real(dp), intent(in) :: A(:, :), Emin, Emax
integer, intent(in) :: M0
real(dp), intent(out) :: lam(:), c(:, :)
integer :: feastparam(64)
real(dp), allocatable :: res(:)
real(dp) :: epsout
integer :: i, loop, info, M, N, LDA
N = size(A, 1)
LDA = N
allocate(res(M0))

call feastinit(feastparam)
feastparam(1)=1 !! change from default value
call dfeast_syev('F',N,A,LDA,feastparam,epsout,loop,Emin,Emax,M0,lam,c, &
    M,res,info)

if (info /= 0) then
    print *,'FEAST OUTPUT INFO', info
    call stop_error("info /= 0")
end if

print *, '# Search interval [Emin,Emax]', Emin, Emax
print *, '# mode found/subspace', M, M0
print *, '# iterations', loop
print *, 'TRACE', sum(lam(1:M))
print *, 'Relative error on the Trace', epsout
print *, 'Eigenvalues/Residuals'
do i = 1, M
    print *, i, lam(i), res(i)
end do

if (M < M0) call stop_error("M < M0")

end subroutine

end module
