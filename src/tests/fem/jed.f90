module jed_functions
use types, only: dp
use constants, only: pi
use linalg, only: eye, solve, eig, inv, eigvals
implicit none
private
public logg

contains

real(dp) function logg(m) result(error)
integer, intent(in) :: m
real(dp) :: D(m+2, m+2), D_(m+2, m+2), xb(m+2), D2(m, m), xd(m), sinx(m), &
    xexact(m, m), b(m, m), x(m,m), lam(m), R(m, m), Diag(m, m), Rinv(m, M)
complex(dp) :: lam_complex(m), R_complex(m, m)
integer :: i, j
call Cheb(m+1, D, xb);
xb = 0.5_dp + 0.5_dp * xb
D_ = -4*matmul(D, D)
D2 = D_(2:size(D,1)-1, 2:size(D,2)-1)
xd = xb(2:size(xb)-1)
sinx(:) = sin(pi*xd)
forall (i=1:m, j=1:m) xexact(i,j) = sinx(i)*sinx(j)
b = 2*pi**2*xexact

call eig(D2, lam_complex, R_complex)
lam = real(lam_complex, dp)
R = real(R_complex, dp)
forall(i=1:m, j=1:m) Diag(i, j) = 1/(lam(i)+lam(j))
Rinv = inv(R)
Diag = Diag * matmul(matmul(Rinv, b), transpose(Rinv))
x = matmul(matmul(R, Diag), transpose(R))

error = norm2(x - xexact)
end function

subroutine Cheb(N, D, x)
integer, intent(in) :: N
real(dp), intent(out) :: D(:, :), x(:)
real(dp) :: c(N+1), X_(N+1, N+1), dX(N+1, N+1)
integer :: i, j
if (N == 0) then; D = 0; x = 1; return; end if
x = cos(pi*[(i,i=0,N)]/N)
c(:) = [2, [(1,i=1,N-1)], 2] * (-1)**[(i,i=0,N)]
X_ = spread(x, 2, N+1)
dX = X_ - transpose(X_)
forall(i=1:N+1) dX(i, i) = dX(i, i) + 1
forall(i=1:N+1, j=1:N+1) D(i, j) = c(i) / c(j)
D = D / dX
forall(i=1:N+1) D(i, i) = D(i, i) - sum(D(i, :))
end subroutine

function kron(A, B) result(C)
real(dp), intent(in) :: A(:, :), B(:, :)
real(dp) :: C(size(A, 1)*size(B, 1), size(A, 2)*size(B, 2))
integer :: i, j
forall(i=1:size(A, 1), j=1:size(A, 2))
    C(size(B,1)*(i-1)+1:size(B,1)*i, size(B,2)*(j-1)+1:size(B,2)*j)=A(i,j)*B
end forall
end function

end module

program jed
use types, only: dp
use jed_functions, only: logg
implicit none
real(dp) :: error
integer :: i, N
real(dp) :: t1, t2

N = 100000
call cpu_time(t1)
do i = 1, N
    error = logg(7)
end do
call cpu_time(t2)
print *, "error = ", error
print "('Total time (us): ', f12.8)", 1e6_dp * (t2-t1) / N
end program
