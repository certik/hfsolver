module fourier

! Calculates the discrete Fourier transform using the simple direct but slow
! O(n^2) algorithm.

use types, only: dp
use constants, only: i_, pi
implicit none
private
public dft, idft

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

end module
