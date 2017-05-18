module mixings

! This module contains SCF mixing algorithms.

use types, only: dp
implicit none
private
public mixing_linear, mixing_linear_adapt

interface
    subroutine R_function(x, y, E)
    ! Computes y = R(x)
    import :: dp
    implicit none
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:), E
    end subroutine
end interface

contains

subroutine mixing_linear(myid, n, R, max_iter, alpha, x)
! Finds "x" so that R(x) = 0
integer, intent(in) :: myid, n, max_iter
real(dp), intent(in) :: alpha
procedure(R_function) :: R
! On input: initial estimate x0; On output: "x" satisfies R(x) = 0
real(dp), intent(inout) :: x(n)
real(dp) :: y(n), E
integer :: i
do i = 1, max_iter
    call R(x, y, E)
    x = x + alpha * y
    if (myid == 0) then
        print *, "ITER:", i, E
    end if
end do
end subroutine

subroutine mixing_linear_adapt(myid, n, R, max_iter, alpha, x)
! Finds "x" so that R(x) = 0, uses x0 as the initial estimate
integer, intent(in) :: myid, n, max_iter
real(dp), intent(in) :: alpha
procedure(R_function) :: R
! On input: initial estimate x0; On output: "x" satisfies R(x) = 0
real(dp), intent(inout) :: x(n)

real(dp), parameter :: alpha_max = 1
real(dp), dimension(n) :: R_m, R_mm1, beta
real(dp) :: E
integer :: i, j
beta = alpha
do i = 1, max_iter
    call R(x, R_m, E)
    if (myid == 0) then
        print *, "ITER:", i, E
    end if
    x = x + beta * R_m
    if (i > 1) then
        do j = 1, n
            if (R_mm1(j) * R_m(j) > 0) then
                beta(j) = beta(j) + alpha
                if (beta(j) > alpha_max) beta(j) = alpha_max
            else
                beta(j) = alpha
            end if
        end do
    end if
    R_mm1 = R_m
end do
end subroutine

end module
