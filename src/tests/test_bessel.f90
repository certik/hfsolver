program test_bessel
use types, only: dp
use constants, only: pi
use special_functions, only: Fm, Inu_formula2
use utils, only: assert, init_random
implicit none
integer :: i, n, k
real(dp) :: r
real(dp), allocatable :: yf(:), yr(:), x(:)
real(dp) :: t1, t2, t3

call init_random()

n = 100000
allocate(x(n), yf(n), yr(n))
do k = 0, 4
    print *, "Testing k =", k
    do i = 1, n
        call random_number(r)
        x(i) = r*40
    end do
    call cpu_time(t1)
    do i = 1, n
        yf(i) = f(k, x(i))
    end do
    call cpu_time(t2)
    do i = 1, n
        yr(i) = Inu_formula2(k, x(i))
    end do
    call cpu_time(t3)
    print *, "abs:", maxval(abs(yf-yr))
    print *, "rel:", maxval(abs(yf-yr) / max(abs(yf), abs(yr)))
    print *, "time f(r):", t2-t1
    print *, "time r(r):", t3-t2
    !print *, "speedup:", (t2-t1) / (t3-t2)
    print *
end do

contains

real(dp) pure function f(k, x) result(r)
integer, intent(in) :: k
real(dp), intent(in) :: x
select case (k)
    case (0)
        r = sinh(x)
    case (1)
        r = -sinh(x)/x + cosh(x)
    case (2)
        r = (3/x**2 + 1)*sinh(x) - 3/x*cosh(x)
    case (3)
        r = -(15/x**3 + 6/x)*sinh(x) + (15/x**2 + 1)*cosh(x)
    case (4)
        r = (105/x**4 + 45/x**2 + 1)*sinh(x) - (105/x**3 + 10/x)*cosh(x)
end select
r = r * sqrt(2/(pi*x)) / exp(x)
end function

end program
