program analyze_md
use types, only: dp
use optimize, only: linregress
use utils, only: assert
implicit none

real(dp), allocatable :: x(:), y(:)
real(dp) :: slope, intercept, r, stderr_slope, stderr_intercept
integer :: i

allocate(x(4), y(4))
x = [5.05, 6.75, 3.21, 2.66]
y = [1.65, 26.5, -5.93, 7.96]
call linregress(x, y, slope, intercept, r, stderr_slope, stderr_intercept)
print *, slope, intercept, r, stderr_slope, stderr_intercept
call assert((slope - 5.3935772582_dp) < 1e-9_dp)
call assert((intercept - (-16.2811279161_dp)) < 1e-9_dp)
call assert((r - 0.7244351257_dp) < 1e-9_dp)
call assert((stderr_slope - 3.6290902263_dp) < 1e-9_dp)
call assert((stderr_intercept - 17.0648974636_dp) < 1e-9_dp)
deallocate(x, y)

allocate(x(10), y(10))
x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
y = [0, 2, 4, 3, 5, 7, 6, 8, 10, 9]
call linregress(x, y, slope, intercept, r, stderr_slope, stderr_intercept)
call assert((slope - 1.0181818182_dp) < 1e-9_dp)
call assert((intercept - 0.8181818182_dp) < 1e-9_dp)
call assert((r - 0.9620913858_dp) < 1e-9_dp)
call assert((stderr_slope - 0.1020452015_dp) < 1e-9_dp)
call assert((stderr_intercept - 0.5447723006_dp) < 1e-9_dp)
deallocate(x, y)

allocate(x(100), y(100))
forall(i=1:100) x(i) = i-1
y = x+modulo(x, 3._dp)
call linregress(x, y, slope, intercept, r, stderr_slope, stderr_intercept)
call assert((slope - 1.0001980198_dp) < 1e-9_dp)
call assert((intercept - 0.9801980198_dp) < 1e-9_dp)
call assert((r - 0.9995984406_dp) < 1e-9_dp)
call assert((stderr_slope - 0.0028641364_dp) < 1e-9_dp)
call assert((stderr_intercept - 0.1641202648_dp) < 1e-9_dp)
deallocate(x, y)

end program
