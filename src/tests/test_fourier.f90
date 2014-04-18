program test_fourier
use types, only: dp
use fourier, only: dft, idft, fft, fft_vectorized, fft_pass, fft_pass_inplace, &
        fft_vectorized_inplace, calculate_factors, ifft_pass, fft2_inplace, &
        fft3_inplace, ifft3_inplace
use utils, only: assert, init_random
use constants, only: i_
implicit none

real(dp), allocatable :: x(:)
complex(dp), allocatable :: xdft(:), x2(:, :), x3(:, :, :), x3d(:, :, :)
real(dp) :: tmp
real(dp) :: t1, t2
integer :: n, i, j, k

! Test 1-16
call test_factors(1, [1])
call test_factors(2, [2])
call test_factors(3, [3])
call test_factors(4, [4])
call test_factors(5, [5])
call test_factors(6, [2, 3])
call test_factors(7, [7])
call test_factors(8, [2, 4])
call test_factors(9, [3, 3])
call test_factors(10, [2, 5])
call test_factors(11, [11])
call test_factors(12, [3, 4])
call test_factors(13, [13])
call test_factors(14, [2, 7])
call test_factors(15, [3, 5])
call test_factors(16, [4, 4])

! Test 2^k
call test_factors(32, [2, 4, 4])
call test_factors(64, [4, 4, 4])
call test_factors(128, [2, 4, 4, 4])
call test_factors(256, [4, 4, 4, 4])
call test_factors(512, [2, 4, 4, 4, 4])
call test_factors(1024, [4, 4, 4, 4, 4])
call test_factors(2048, [2, 4, 4, 4, 4, 4])
call test_factors(1024**2/2, [2, 4, 4, 4, 4, 4, 4, 4, 4, 4])
call test_factors(1024**2, [4, 4, 4, 4, 4, 4, 4, 4, 4, 4])

! Test various special values:
call test_factors(24, [2, 3, 4])
call test_factors(100, [4, 5, 5])
call test_factors(121, [11, 11])
call test_factors(169, [13, 13])
call test_factors(289, [17, 17])
call test_factors(343, [7, 7, 7])
call test_factors(384, [2, 3, 4, 4, 4])
call test_factors(1000, [2, 4, 5, 5, 5])
call test_factors(1309, [7, 11, 17])
call test_factors(2401, [7, 7, 7, 7])
call test_factors(3030, [2, 3, 5, 101])
call test_factors(3360, [2, 3, 4, 4, 5, 7])


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

! Test 1-100
do i = 1, 100
    call test_fft_pass(i)
end do
! Test special cases
call test_fft_pass(121)
call test_fft_pass(169)
call test_fft_pass(289, eps=1e-8_dp)
call test_fft_pass(343, eps=1e-8_dp)

allocate(x2(3, 5))
forall(i=1:size(x2, 1), j=1:size(x2, 2)) x2(i, j) = i*j+3*i**2
call fft2_inplace(x2)
call assert(abs(sum(x2)-60) < 1e-9_dp)

allocate(x3(3, 5, 7))
forall(i=1:size(x3, 1), j=1:size(x3, 2), k=1:size(x3, 3))
    x3(i, j, k) = i*j+3*i**2+k*j+k
end forall
call fft3_inplace(x3)
call assert(abs(sum(x3)-630) < 1e-9_dp)
deallocate(x3)

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

call cpu_time(t1)
xdft = ifft_pass(xdft)
call cpu_time(t2)
print *, "ifft_pass"
print *, "time:", (t2-t1)*1000, "ms"
deallocate(x, xdft)

print *, "fft3:"
n = 16
allocate(x3(n, n, n), x3d(n, n, n))
do i = 1, size(x3, 1)
do j = 1, size(x3, 2)
do k = 1, size(x3, 3)
    call random_number(tmp)
    x3(i, j, k) = tmp
end do
end do
end do
x3d = x3
call cpu_time(t1)
call fft3_inplace(x3d)
call cpu_time(t2)
print *, "time:", (t2-t1)*1000, "ms"
call cpu_time(t1)
call ifft3_inplace(x3d)
call cpu_time(t2)
print *, "time:", (t2-t1)*1000, "ms"
call assert(all(abs(x3 - x3d/n**3) < 1e-15_dp))
deallocate(x3, x3d)

contains

subroutine test_factors(n, correct_fac)
integer, intent(in) :: n, correct_fac(:)
integer, allocatable :: fac(:)
call calculate_factors(n, fac)
call assert(all(fac == correct_fac))
end subroutine

subroutine test_fft_pass(n, eps)
integer, intent(in) :: n
real(dp), optional, intent(in) :: eps
real(dp), allocatable :: x(:)
real(dp) :: eps_
if (present(eps)) then
    eps_ = eps
else
    eps_ = 1e-9_dp
end if
allocate(x(n))
forall(i = 1:n) x(i) = i
call assert(all(abs(dft(x) - fft_pass(x)) < eps_))
call assert(all(abs(x - ifft_pass(fft_pass(x))/n) < eps_))
forall(i = 1:n) x(i) = 1._dp / i
call assert(all(abs(dft(x) - fft_pass(x)) < eps_))
call random_number(x)
call assert(all(abs(dft(x) - fft_pass(x)) < eps_))
end subroutine

end program
