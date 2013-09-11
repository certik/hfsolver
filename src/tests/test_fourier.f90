program test_fourier
use types, only: dp
use fourier, only: dft, idft, fft, fft_vectorized, fft_pass, fft_pass_inplace
use utils, only: assert, init_random
use constants, only: i_
implicit none

real(dp), allocatable :: x(:)
complex(dp), allocatable :: xdft(:)
real(dp) :: t1, t2
integer :: n, i

call assert(all(abs(dft([1._dp]) - [1]) < 1e-12_dp))
call assert(all(abs(idft([1._dp+0*i_]) - [1]) < 1e-12_dp))
call assert(all(abs(dft([1._dp, 2._dp]) - [3, -1]) < 1e-12_dp))
call assert(all(abs(idft([3._dp+0*i_, -1._dp + 0*i_]) - [1, 2]) < 1e-12_dp))

allocate(x(4), xdft(4))
x = [1, 2, 3, 4]
xdft = [10 + 0*i_, -2+2*i_, -2+0*i_, -2-2*i_]
call assert(all(abs(dft(x) - xdft) < 1e-12_dp))
call assert(all(abs(fft(x) - xdft) < 1e-12_dp))
call assert(all(abs(fft_vectorized(x) - xdft) < 1e-12_dp))
call assert(all(abs(fft_pass(x) - xdft) < 1e-12_dp))
call assert(all(abs(idft(xdft) - x) < 1e-12_dp))
deallocate(x, xdft)

allocate(x(5))
x = [1, 2, 3, 4, 5]
call assert(all(abs(idft(dft(x)) - x) < 1e-12_dp))
deallocate(x)

n = 8
allocate(x(n))
forall(i = 1:n) x(i) = i
call assert(all(abs(idft(dft(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft_vectorized(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft_pass(x)) - x) < 1e-10_dp))
deallocate(x)

n = 16
allocate(x(n))
forall(i = 1:n) x(i) = i
call assert(all(abs(idft(dft(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft_vectorized(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft_pass(x)) - x) < 1e-10_dp))
deallocate(x)

n = 32
allocate(x(n))
forall(i = 1:n) x(i) = i
call assert(all(abs(idft(dft(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft_vectorized(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft_pass(x)) - x) < 1e-10_dp))
deallocate(x)

n = 64
allocate(x(n))
forall(i = 1:n) x(i) = i
call assert(all(abs(idft(dft(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft_vectorized(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft_pass(x)) - x) < 1e-10_dp))
deallocate(x)

n = 128
allocate(x(n))
forall(i = 1:n) x(i) = i
call assert(all(abs(idft(dft(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft_vectorized(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft_pass(x)) - x) < 1e-10_dp))
deallocate(x)

n = 256
allocate(x(n))
forall(i = 1:n) x(i) = i
call assert(all(abs(idft(dft(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft_vectorized(x)) - x) < 1e-10_dp))
call assert(all(abs(idft(fft_pass(x)) - x) < 1e-10_dp))
deallocate(x)

n = 400
allocate(x(n))
forall(i = 1:n) x(i) = i
call assert(all(abs(idft(dft(x)) - x) < 1e-10_dp))
deallocate(x)

n = 1024
call init_random()
allocate(x(n), xdft(n))
call random_number(x)
call cpu_time(t1)
xdft = fft(x)
call cpu_time(t2)
print *, "fft"
print *, "time:", (t2-t1)*1000, "ms"

call cpu_time(t1)
xdft = fft_vectorized(x)
call cpu_time(t2)
print *, "fft_vectorized"
print *, "time:", (t2-t1)*1000, "ms"

call cpu_time(t1)
xdft = fft_pass(x)
call cpu_time(t2)
print *, "fft_pass"
print *, "time:", (t2-t1)*1000, "ms"
deallocate(x, xdft)

end program
