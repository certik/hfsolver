module isolve
use types, only: dp
use utils, only: stop_error
implicit none
private
public solve_cg

contains

function solve_cg(A, b, x0, tol, maxiter) result(x)
! solves a system of equations A x = b with one right hand side
real(dp), intent(in) :: A(:, :)  ! coefficient matrix A
real(dp), intent(in) :: b(:)  ! right-hand-side A x = b
real(dp), intent(in) :: x0(:) ! Initial guess
real(dp), intent(in) :: tol ! Tolerance in residual
integer, intent(in) :: maxiter ! Maximum number of iterations
real(dp) :: x(size(b)) ! solution
real(dp), dimension(size(b)) :: r, p, Ap
real(dp) :: r2, r2old, alpha
integer :: i

x = x0
r = b - matmul(A, x)
p = r
r2old = dot_product(r, r)    ! ||A x - b||^2
print *, "Conjugate Gradient solver"
print *, "Iter    Residual ||A x - b||"
do i = 1, maxiter
    Ap = matmul(A, p)
    alpha = r2old / dot_product(p, Ap)
    x = x + alpha * p
    r = r - alpha * Ap
    r2 = dot_product(r, r)  ! Good approximation for ||A x - b||^2

    print "(i4, '      ', es10.2)", i, sqrt(r2)
    if (sqrt(r2) < tol) then
        r = matmul(A, x) - b
        write(*, '(1x,a,es10.2)') "Solution vector residual ||A x - b||/||bv||: ",  sqrt(dot_product(r, r) / dot_product(b, b))
        return
    end if

    p = r + r2 / r2old * p
    r2old = r2
end do
call stop_error("Solution did not converge.")
end function

end module
