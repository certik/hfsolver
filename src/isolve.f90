module isolve
use types, only: dp
use utils, only: stop_error
use sparse, only: csr_matvec
implicit none
private
public solve_cg

contains

function solve_cg(Ap, Aj, Ax, b, x0, tol, maxiter) result(x)
! solves a system of equations A x = b with one right hand side
integer, intent(in) :: Ap(:), Aj(:)  ! coefficient matrix A in CSR format
real(dp), intent(in) :: Ax(:)        ! coefficient matrix A in CSR format
real(dp), intent(in) :: b(:)  ! right-hand-side A x = b
real(dp), intent(in) :: x0(:) ! Initial guess
real(dp), intent(in) :: tol ! Tolerance in residual
integer, intent(in) :: maxiter ! Maximum number of iterations
real(dp) :: x(size(b)) ! solution
real(dp), dimension(size(b)) :: r, p, Ap_, z
real(dp) :: r2, r2old, alpha, res_norm
integer :: i

x = x0
r = b - csr_matvec(Ap, Aj, Ax, x)
z = precond(r)
p = z
r2old = dot_product(r, z)
print *, "Conjugate Gradient solver"
print *, "Iter    Residual ||A x - b||"
do i = 1, maxiter
    Ap_ = csr_matvec(Ap, Aj, Ax, x)
    alpha = r2old / dot_product(p, Ap_)
    x = x + alpha * p
    r = r - alpha * Ap_
    res_norm = sqrt(dot_product(r, r))  ! Good approximation to ||A x - b||

    print "(i4, '      ', es10.2)", i, res_norm
    if (res_norm < tol) then
        r = csr_matvec(Ap, Aj, Ax, x) - b
        write(*, '(1x,a,es10.2)') "Solution vector residual ||A x - b||/||bv||: ",  sqrt(dot_product(r, r) / dot_product(b, b))
        return
    end if

    z = precond(r)
    r2 = dot_product(r, z)
    p = z + r2 / r2old * p
    r2old = r2
end do
call stop_error("Solution did not converge.")

contains

    function precond(x) result(y)
    ! Calculates y = M^-1 x
    real(dp), intent(in) :: x(:)
    real(dp) :: y(size(x))
    ! TODO: allow to retrieve elements (i, j) from the CSR matrix (write a
    ! subroutine and reenble this:
    !integer :: i
    ! No preconditioning:
    y = x
    ! Jacobi normalization: M = diag(A):
    !do i = 1, size(x)
    !    y(i) = x(i) / A(i, i)
    !end do
    end function

end function

end module
