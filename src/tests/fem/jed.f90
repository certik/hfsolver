module jed_functions
use types, only: dp
use constants, only: pi
use linalg, only: eye, solve, eig, inv, eigvals
implicit none
private
public logg

contains

real(dp) function logg(m, solvetype) result(error)
integer, intent(in) :: m, solvetype
real(dp) :: D(m+2, m+2), xb(m+2), D2(m, m), xd(m), sinx(m, 1), xexact(m**2), &
    b(m**2), L(m**2, m**2), x(m**2), lam(m), R(m, m), Diag(m, m), Rinv(m, M)
complex(dp) :: lam_complex(m), R_complex(m, m)
integer :: i, j
call Cheb(m+1, D, xb);
D = 2*D
xb = 0.5_dp + 0.5_dp * xb
D = -matmul(D, D)
D2 = D(2:size(D,1)-1, 2:size(D,2)-1)
xd = xb(2:size(xb)-1)
sinx(:, 1) = sin(pi*xd)
xexact = reshape(matmul(sinx, transpose(sinx)), [size(xexact)])
b = 2*pi**2*xexact
select case (solvetype)
    case (1)
        L = kron(D2, eye(m)) + kron(eye(m), D2)
        x = solve(L, b)
    case (2)
        call eig(D2, lam_complex, R_complex)
        lam = real(lam_complex, dp)
        R = real(R_complex, dp)
        forall(i=1:m, j=1:m) Diag(i, j) = 1/(lam(i)+lam(j))
        Rinv = inv(R)
        Diag = Diag * matmul(matmul(Rinv, reshape(b, [m, m])), transpose(Rinv))
        x = reshape(matmul(matmul(R, Diag), transpose(R)), [size(x)])
end select
error = norm2(x - xexact)
end function

subroutine Cheb(N, D, x)
integer, intent(in) :: N
real(dp), intent(out) :: D(:, :), x(:)
real(dp) :: c(N+1, 1), X_(N+1, N+1), dX(N+1, N+1)
integer :: i
if (N == 0) then; D = 0; x = 1; return; end if
x = cos(pi*[(i,i=0,N)]/N)
c(:, 1) = [2, [(1,i=1,N-1)], 2] * (-1)**[(i,i=0,N)]
X_ = spread(x, 2, N+1)
dX = X_ - transpose(X_)
D = (matmul(c, transpose(1/c))) / (dX + eye(N+1))
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
    error = logg(7, 2)
end do
call cpu_time(t2)
print *, "error = ", error
print "('Total time: ', f10.8)", (t2-t1) / N
end program
