program test_fftw
use iso_c_binding, only: c_ptr
use types, only: dp
use constants, only: i_
use fftw, only: fftw_plan_dft_1d, fftw_plan_dft_3d, FFTW_FORWARD, &
    FFTW_ESTIMATE, fftw_execute_dft, fftw_destroy_plan, &
    FFTW_MEASURE, alloc1d, alloc3d, free
implicit none

type(c_ptr) :: plan
complex(dp), dimension(4, 8, 16) :: y3, ydft3
integer :: n
complex(dp), pointer :: x(:), xdft(:), x3(:, :, :), xdft3(:, :, :)
real(dp) :: t1, t2, t3

! This works, but is slow due to y,ydft not being aligned to exploit SIMD
! The dimensions (4, 8, 16) must be reversed:
plan = fftw_plan_dft_3d(16, 8, 4, y3, ydft3, FFTW_FORWARD, FFTW_ESTIMATE)
y3 = 0
call fftw_execute_dft(plan, y3, ydft3)
call fftw_destroy_plan(plan)

n = 1024**2

! This is fast
print *, "1D FFT of size n=", n, "with FFTW allocation"
x => alloc1d(n)
xdft => alloc1d(n)
x = 0
call cpu_time(t1)
!plan = fftw_plan_dft_1d(n, x, xdft, FFTW_FORWARD, FFTW_ESTIMATE)
plan = fftw_plan_dft_1d(n, x, xdft, FFTW_FORWARD, FFTW_MEASURE)
call cpu_time(t2)
call fftw_execute_dft(plan, x, xdft)
call cpu_time(t3)
print *, "Total time:", (t3-t1)*1000, "ms"
print *, "init:      ", (t2-t1)*1000, "ms"
print *, "calc:      ", (t3-t2)*1000, "ms"
call fftw_destroy_plan(plan)
call free(x)
call free(xdft)

n = 256

! This is fast
print *, "1D FFT of size n=", n, "^3  with FFTW allocation"
x3 => alloc3d(n, n, n)
xdft3 => alloc3d(n, n, n)
x3 = 0
call cpu_time(t1)
plan = fftw_plan_dft_3d(n, n, n, x3, xdft3, FFTW_FORWARD, FFTW_MEASURE)
call cpu_time(t2)
call fftw_execute_dft(plan, x3, xdft3)
call cpu_time(t3)
print *, "Total time:", (t3-t1)*1000, "ms"
print *, "init:      ", (t2-t1)*1000, "ms"
print *, "calc:      ", (t3-t2)*1000, "ms"
call fftw_destroy_plan(plan)
call free(x3)
call free(xdft3)
end program
