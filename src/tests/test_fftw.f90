program test_fftw
use iso_c_binding, only: c_ptr, c_f_pointer, c_size_t, c_loc
use types, only: dp
use constants, only: i_
use fftw, only: fftw_plan_dft_1d, fftw_plan_dft_3d, FFTW_FORWARD, &
    FFTW_ESTIMATE, fftw_execute_dft, fftw_destroy_plan, fftw_alloc_complex, &
    fftw_free, FFTW_MEASURE
implicit none

type(c_ptr) :: plan
complex(dp), dimension(4, 8, 16) :: y3, ydft3
complex(dp), allocatable :: x(:), xdft(:), x3(:, :, :), xdft3(:, :, :)
integer :: n
complex(dp), pointer :: px(:), pxdft(:), px3(:, :, :), pxdft3(:, :, :)
type(c_ptr) :: p1, p2
real(dp) :: t1, t2, t3

! This works, but is slow due to y,ydft not being aligned to exploit SIMD
! The dimensions (4, 8, 16) must be reversed:
plan = fftw_plan_dft_3d(16, 8, 4, y3, ydft3, FFTW_FORWARD, FFTW_ESTIMATE)
call fftw_execute_dft(plan, y3, ydft3)
call fftw_destroy_plan(plan)

n = 1024
! This works, but is slow due to x,xdft not being aligned to exploit SIMD
print *, "1D FFT of size n=", n, "with Fortran allocation"
allocate(x(n), xdft(n))
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
deallocate(x, xdft)

! This is fast
print *, "1D FFT of size n=", n, "with FFTW allocation"
p1 = fftw_alloc_complex(int(n, c_size_t))
call c_f_pointer(p1, px, [n])
p2 = fftw_alloc_complex(int(n, c_size_t))
call c_f_pointer(p2, pxdft, [n])
call cpu_time(t1)
!plan = fftw_plan_dft_1d(n, px, pxdft, FFTW_FORWARD, FFTW_ESTIMATE)
plan = fftw_plan_dft_1d(n, px, pxdft, FFTW_FORWARD, FFTW_MEASURE)
call cpu_time(t2)
call fftw_execute_dft(plan, px, pxdft)
call cpu_time(t3)
print *, "Total time:", (t3-t1)*1000, "ms"
print *, "init:      ", (t2-t1)*1000, "ms"
print *, "calc:      ", (t3-t2)*1000, "ms"
call fftw_destroy_plan(plan)
call fftw_free(p1)
call fftw_free(p2)

n = 16
! This works, but is slow due to x,xdft not being aligned to exploit SIMD
print *, "3D FFT of size n=", n, "^3  with Fortran allocation"
allocate(x3(n, n, n), xdft3(n, n, n))
call cpu_time(t1)
plan = fftw_plan_dft_3d(n, n, n, x3, xdft3, FFTW_FORWARD, FFTW_MEASURE)
call cpu_time(t2)
call fftw_execute_dft(plan, x3, xdft3)
call cpu_time(t3)
print *, "Total time:", (t3-t1)*1000, "ms"
print *, "init:      ", (t2-t1)*1000, "ms"
print *, "calc:      ", (t3-t2)*1000, "ms"
call fftw_destroy_plan(plan)
deallocate(x3, xdft3)

! This is fast
print *, "1D FFT of size n=", n, "^3  with FFTW allocation"
p1 = fftw_alloc_complex(int(n**3, c_size_t))
call c_f_pointer(p1, px3, [n, n, n])
p2 = fftw_alloc_complex(int(n**3, c_size_t))
call c_f_pointer(p2, pxdft3, [n, n, n])
call cpu_time(t1)
plan = fftw_plan_dft_3d(n, n, n, px3, pxdft3, FFTW_FORWARD, FFTW_MEASURE)
call cpu_time(t2)
call fftw_execute_dft(plan, px3, pxdft3)
call cpu_time(t3)
print *, "Total time:", (t3-t1)*1000, "ms"
print *, "init:      ", (t2-t1)*1000, "ms"
print *, "calc:      ", (t3-t2)*1000, "ms"
call fftw_destroy_plan(plan)
call fftw_free(p1)
call fftw_free(p2)
end program
