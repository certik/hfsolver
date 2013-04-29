module isolve
use types, only: dp
implicit none
private
public solve_cg

contains

function solve_cg(A, b) result(x)
! solves a system of equations A x = b with one right hand side
real(dp), intent(in) :: A(:,:)  ! coefficient matrix A
real(dp), intent(in) :: b(:)  ! right-hand-side A x = b
real(dp) :: x(size(b))
real(dp), dimension(size(b)) :: r, p, q
real(dp) :: rho, rho_prev, alpha, beta
integer :: i, max_iter

x = 0
r = b
rho = dot_product(r, r)
p = r
max_iter = 100
do i = 1, max_iter
    q = matmul(A, p)
    alpha = rho / dot_product(p, q)
    x = x + alpha * p
    r = r - alpha * q
    rho_prev = rho
    rho = dot_product(r, r)
    beta = rho / rho_prev
    p = r + beta * p
    print *, i, rho
end do
end function

end module
