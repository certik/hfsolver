module random

use types, only: dp
implicit none
private
public randn

interface randn
    module procedure randn_scalar
    module procedure randn_vector
    module procedure randn_matrix
    module procedure randn_vector_n
end interface

contains

subroutine randn_scalar(x)
! Returns a psuedorandom scalar drawn from the standard normal distribution.
!
! [1] Marsaglia, G., & Bray, T. A. (1964). A Convenient Method for Generating
!       Normal Variables. SIAM Review, 6(3), 260–264.
real(dp), intent(out) :: x
logical, save :: first = .true.
real(dp), save :: u(2)
real(dp) :: r2
if (first) then
    do
        call random_number(u)
        u = 2*u-1
        r2 = sum(u**2)
        if (r2 < 1 .and. r2 > 0) exit
    end do
    u = u * sqrt(-2*log(r2)/r2)
    x = u(1)
else
    x = u(2)
end if
first = .not. first
end subroutine

subroutine randn_vector_n(n, x)
integer, intent(in) :: n
real(dp), intent(out) :: x(n)
integer :: i
do i = 1, size(x)
    call randn(x(i))
end do
end subroutine

subroutine randn_vector(x)
real(dp), intent(out) :: x(:)
call randn_vector_n(size(x), x)
end subroutine

subroutine randn_matrix(x)
real(dp), intent(out) :: x(:, :)
call randn_vector_n(size(x), x)
end subroutine

end module
