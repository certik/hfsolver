module fourier

! Calculates the discrete Fourier transform using the simple direct but slow
! O(n^2) algorithm.

use types, only: dp
use constants, only: i_, pi
use utils, only: stop_error, assert
implicit none
private
public dft, idft, fft, fft_vectorized, fft_pass, fft_pass_inplace

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

subroutine dft_vec(Ns, N, x, p)
! Compute the one-dimensional discrete Fourier transform on each row of 'x'
! separately.
integer, intent(in) :: Ns, N
real(dp), intent(in) :: x(Ns, N)
complex(dp), intent(out) :: p(Ns, N)
complex(dp) :: F(N, N)
integer :: i, j
forall(i=0:N-1, j=0:N-1, i >= j) F(i+1, j+1) = exp(-2*pi*i_*i*j / N)
forall(i=1:N, j=1:N, i < j) F(i, j) = F(j, i)
p = matmul(x, F)
end subroutine

subroutine fft_step(Ns, Nmin, x, p)
integer, intent(in) :: Ns, Nmin
complex(dp), intent(in) :: x(Ns, Nmin)
complex(dp), intent(out) :: p(Ns/2, Nmin*2)
complex(dp) :: tmp
integer :: i
do i = 1, Nmin
    !tmp = exp(-pi*i_*(i-1)/Nmin)
    ! The same as the previous line, just faster:
    tmp = cos(pi*(i-1)/Nmin) - i_*sin(pi*(i-1)/Nmin)
    p(:,      i) = x(:Ns/2, i) + tmp * x(Ns/2+1:, i)
    p(:, Nmin+i) = x(:Ns/2, i) - tmp * x(Ns/2+1:, i)
end do
end subroutine

function fft_vectorized(x) result(p)
! A vectorized, non-recursive version of the Cooley-Tukey FFT
real(dp), intent(in), target :: x(:)
complex(dp), target :: p(size(x))
complex(dp), target :: tmp(size(x))
integer :: N, Nmin, Ns
logical :: p_is_result
N = size(x)
if (iand(N, N-1) /= 0) call stop_error("size of x must be a power of 2")
Nmin = min(N, 4)
Ns = N / Nmin
call dft_vec(Ns, Nmin, x, p)
p_is_result = .true.
do while (Nmin < N)
    if (p_is_result) then
        call fft_step(Ns, Nmin, p, tmp)
    else
        call fft_step(Ns, Nmin, tmp, p)
    end if
    Nmin = Nmin * 2
    Ns = Ns / 2
    p_is_result = .not. p_is_result
end do
if (.not. p_is_result) p = tmp
end function

subroutine precalculate_coeffs(wa)
! Precalculates all cos/sin factors
complex(dp), intent(out) :: wa(:)
integer :: n, k, i, idx
n = size(wa)
k = n / 2
idx = 1
do while (k > 0)
    wa(idx) = 1
    do i = 1, k
        idx = idx + 1
        ! Equivalent to exp(-i_*i*pi/k) but faster:
        wa(idx) = cos(i*pi/k) - i_ * sin(i*pi/k)
    end do
    k = k/2
end do
end subroutine

subroutine passf2(IDO, L1, CC, CH, WA1)
! FFT pass of factor 2
integer, intent(in) :: IDO, L1
complex(dp), intent(in) :: CC(IDO, 2, L1), WA1(:)
complex(dp), intent(out) :: CH(IDO, L1, 2)
integer :: I, K
do K = 1, L1
    do I = 1, IDO
        CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
        CH(I,K,2) = WA1(I)*(CC(I,1,K)-CC(I,2,K))
    end do
end do
end subroutine

subroutine fft_pass_inplace(x)
complex(dp), intent(inout) :: x(:)
complex(dp), dimension(size(x)) :: angles, CH
integer :: n, nf
integer :: IDOT, IW, K1, L1
logical :: CH_is_result
n = size(x)
if (iand(n, n-1) /= 0) call stop_error("size of x must be a power of 2")
nf = int(log(real(n, dp)) / log(2._dp) + 0.5_dp)
! The above is robust, but let's check it just in case:
call assert(2**nf == n)
call precalculate_coeffs(angles)
CH_is_result = .true.
L1 = 1
IW = 1
do K1 = 1, NF
    IDOT = N/L1/2
    if (CH_is_result) then
        CALL passf2(IDOT,L1,x,CH,angles(IW:IW+IDOT))
    else
        CALL passf2(IDOT,L1,CH,x,angles(IW:IW+IDOT))
    end if
    CH_is_result = .not. CH_is_result
    L1 = 2*L1
    IW = IW+IDOT
end do
if (.not. CH_is_result) then
    x = CH
end if
end subroutine

function fft_pass(x) result(p)
real(dp), intent(in) :: x(:)
complex(dp) :: p(size(x))
p = x
call fft_pass_inplace(p)
end function

end module
