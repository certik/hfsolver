module isolve
use types, only: dp
use utils, only: stop_error
use sparse, only: csr_matvec, csr_getvalue
implicit none
private
public solve_cg

contains

function solve_cg(Ap, Aj, Ax, b, x0, tol, maxiter, verbose) result(x)
! solves a system of equations A x = b with one right hand side
integer, intent(in) :: Ap(:), Aj(:)  ! coefficient matrix A in CSR format
real(dp), intent(in) :: Ax(:)        ! coefficient matrix A in CSR format
real(dp), intent(in) :: b(:)  ! right-hand-side A x = b
real(dp), intent(in) :: x0(:) ! Initial guess
real(dp), intent(in) :: tol ! Tolerance in residual
integer, intent(in) :: maxiter ! Maximum number of iterations
logical, intent(in), optional :: verbose
real(dp) :: x(size(b)) ! solution
real(dp), dimension(size(b)) :: r, p, Ap_, z
real(dp) :: r2, r2old, alpha, res_norm
real(dp) :: A_diag(size(x0))
integer :: i
logical :: verbose_
verbose_ = .false.
if (present(verbose)) verbose_ = verbose
call precond_init()

x = x0
r = b - csr_matvec(Ap, Aj, Ax, x)
z = precond(r)
p = z
r2old = dot_product(r, z)
if (verbose_) then
    print *, "Conjugate Gradient solver"
    print *, "Iter    Residual ||A x - b||"
end if
do i = 1, maxiter
    Ap_ = csr_matvec(Ap, Aj, Ax, p)
    alpha = r2old / dot_product(p, Ap_)
    x = x + alpha * p
    r = r - alpha * Ap_
    res_norm = sqrt(dot_product(r, r))  ! Good approximation to ||A x - b||

    if (verbose_) then
        print "(i4, '      ', es10.2)", i, res_norm
    end if
    if (res_norm < tol) then
        if (verbose_) then
            r = csr_matvec(Ap, Aj, Ax, x) - b
            write(*, '(1x,a,es10.2)') "Solution vector residual ||A x - b||/||bv||: ",  sqrt(dot_product(r, r) / dot_product(b, b))
        end if
        return
    end if

    z = precond(r)
    r2 = dot_product(r, z)
    p = z + r2 / r2old * p
    r2old = r2
end do
call stop_error("Solution did not converge.")

contains

    subroutine precond_init()
    ! Precalculate the diagonal of the A matrix
    integer :: i
    do i = 1, size(x)
        A_diag(i) = csr_getvalue(Ap, Aj, Ax, i, i)
    end do
    end subroutine

    function precond(x) result(y)
    ! Calculates y = M^-1 x
    real(dp), intent(in) :: x(:)
    real(dp) :: y(size(x))
    ! No preconditioning:
    !y = x
    ! Jacobi preconditioning: M = diag(A):
    y = x / A_diag
    end function

end function

end module
