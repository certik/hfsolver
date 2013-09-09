module fourier

! Calculates the discrete Fourier transform using the simple direct but slow
! O(n^2) algorithm.

use types, only: dp
use constants, only: i_, pi
use utils, only: stop_error
implicit none
private
public dft, idft, fft, fft_vectorized

contains

function dft(x) result(p)
! Compute the one-dimensional discrete Fourier transform
real(dp), intent(in) :: x(:)
complex(dp) :: p(size(x))
complex(dp) :: F(size(x), size(x))
integer :: N, i, j
N = size(x)
forall(i=0:N-1, j=0:N-1, i >= j) F(i+1, j+1) = exp(-2*pi*i_*i*j / N)
forall(i=1:N, j=1:N, i < j) F(i, j) = F(j, i)
p = matmul(F, x)
end function

function idft(p) result(x)
! Compute the one-dimensional inverse discrete Fourier transform
! The normalization is such that idft(dft(x)) == x to within numerical
! accuracy.
complex(dp), intent(in) :: p(:)
complex(dp) :: x(size(p))
complex(dp) :: F(size(p), size(p))
integer :: N, i, j
N = size(p)
forall(i=0:N-1, j=0:N-1, i >= j) F(i+1, j+1) = exp(+2*pi*i_*i*j/ N)
forall(i=1:N, j=1:N, i < j) F(i, j) = F(j, i)
x = matmul(F, p) / N
end function

recursive function fft(x) result(p)
! A recursive implementation of the 1D Cooley-Tukey FFT
real(dp), intent(in) :: x(:)
complex(dp) :: p(size(x)), X_even(size(x)/2), X_odd(size(x)/2)
complex(dp) :: factor(size(x))
integer :: N, i
N = size(x)
if (iand(N, N-1) /= 0) call stop_error("size of x must be a power of 2")
if (N <= 4) then
    p = dft(x)
else
    X_even = fft(x(::2))
    X_odd = fft(x(2::2))
    forall(i=0:N-1) factor(i+1) = exp(-2*pi*i_*i/N)
    p = [X_even + factor(:N / 2) * X_odd, X_even + factor(N / 2+1:) * X_odd]
end if
end function

function dft_vec(x) result(p)
! Compute the one-dimensional discrete Fourier transform on each column of 'x'
! separately.
real(dp), intent(in) :: x(:, :)
complex(dp) :: p(size(x, 1), size(x, 2))
complex(dp) :: F(size(x, 1), size(x, 1))
integer :: N, i, j
N = size(x, 1)
forall(i=0:N-1, j=0:N-1, i >= j) F(i+1, j+1) = exp(-2*pi*i_*i*j / N)
forall(i=1:N, j=1:N, i < j) F(i, j) = F(j, i)
p = matmul(F, x)
end function

subroutine fft_step(x, p)
complex(dp), intent(in) :: x(:, :)
complex(dp), intent(out) :: p(:, :)
complex(dp), target :: fac(size(x, 1))
integer :: Nmin, Ns, i
Nmin = size(x, 1)
Ns = size(x, 2)
forall(i=0:Nmin-1) fac(i+1) = exp(-pi*i_*i/Nmin)
forall(i=1:Ns/2) p(:Nmin,   i) = x(:, i) + fac * x(:, Ns/2+i)
forall(i=1:Ns/2) p(Nmin+1:, i) = x(:, i) - fac * x(:, Ns/2+i)
end subroutine

function fft_vectorized(x) result(p)
! A vectorized, non-recursive version of the Cooley-Tukey FFT
real(dp), intent(in) :: x(:)
complex(dp), target :: p(size(x))
complex(dp), target :: tmp(size(x))
complex(dp), pointer :: p1(:, :), p2(:, :)
integer :: N, Nmin, Ns
logical :: p_is_result
N = size(x)
if (iand(N, N-1) /= 0) call stop_error("size of x must be a power of 2")
Nmin = min(N, 32)
Ns = N / Nmin
p = reshape(dft_vec(reshape(x, [Nmin, Ns], order=[2, 1])), [N])
p_is_result = .true.
do while (Nmin < N)
    if (p_is_result) then
        p1(1:Nmin, 1:Ns) => p
    else
        p1(1:Nmin, 1:Ns) => tmp
    end if
    Nmin = Nmin * 2
    Ns = Ns / 2
    if (p_is_result) then
        p2(1:Nmin, 1:Ns) => tmp
    else
        p2(1:Nmin, 1:Ns) => p
    end if
    call fft_step(p1, p2)
    p_is_result = .not. p_is_result
end do
if (.not. p_is_result) p = tmp
end function

end module
