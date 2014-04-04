program test_random
use types, only: dp
use random, only: randn
use utils, only: assert, init_random
implicit none

real(dp) :: a, x(100000), m(10, 10)
call init_random()

! Uniform distribution [0, 1), average 0.5, variance 1/12
call random_number(a) ! Test scalar
call random_number(x) ! Test vector
call random_number(m) ! Test matrix
call tests(x, 0.5_dp, 1._dp/12)

! Standard normal distribution, average 0, variance 1:
call randn(a) ! Test scalar
call randn(x) ! Test vector
call randn(m) ! Test matrix
call tests(x, 0._dp, 1._dp)

contains

    subroutine tests(x, avg0, var0)
    ! We will test for 5 sigma (this test should pass in 99.9999426697% cases,
    ! i.e. fail in one part of 1,744,278):
    real(dp), intent(in) :: x(:), avg0, var0
    real(dp) :: avg, sigma_avg, var, sigma_var
    integer :: n
    n = size(x)

    ! Mean test
    avg = sum(x)/n
    sigma_avg = 1/sqrt(real(n, dp))
    print *, "avg =", avg, "sigma =", sigma_avg
    call assert(abs(avg-avg0) < 5*sigma_avg)

    ! Variance test
    var = sum((x-avg)**2)/n
    sigma_var = sqrt(2._dp/(n-1))
    print *, "var =", var, "sigma =", sigma_var
    call assert(abs(var-var0) < 5*sigma_var)
    end subroutine

end program
