module pofdft_fft
use types, only: dp
use ffte, only: dp_ffte
use pffte, only: pfft3_init, pfft3, pifft3
implicit none
private
public pfft3_init, preal2fourier, pfourier2real


contains

subroutine preal2fourier(x, xG, commy, commz, Ng, nsub)
! Calculates Discrete Fourier Transform, with the same normalization as the
! Fourier Transform for periodic systems (which is Fourier Series).
! Parallel version.
complex(dp), intent(in) :: x(:, :, :)
complex(dp), intent(out) :: xG(:, :, :)
integer, intent(in) :: commy, commz ! communicators in y, z directions
integer, intent(in) :: Ng(:) ! Total (global) number of PW in each direction
integer, intent(in) :: nsub(:) ! Number of subdomains in each direction
! Temporary input array, will get trashed by pfft3
complex(dp) :: tmp(size(x,1), size(x,2), size(x,3))
tmp = x
! Calculates sum_{n=0}^{N-1} e^{-2*pi*i*k*n/N}*x(n)
call pfft3(tmp, xG, commy, commz, Ng, nsub)
xG = xG / product(Ng)     ! The proper normalization is to divide by N
end subroutine

subroutine pfourier2real(xG, x, commy, commz, Ng, nsub)
! Calculates Inverse Discrete Fourier Transform. xG must follow the same
! normalization as defined by real2fourier(), i.e. in the following calls, 'x2'
! will be equal to 'x':
! call real2fourier(x, xG)
! call fourier2real(xG, x2)
! Parallel version.
complex(dp), intent(in) :: xG(:, :, :)
complex(dp), intent(out) :: x(:, :, :)
integer, intent(in) :: commy, commz ! communicators in y, z directions
integer, intent(in) :: Ng(:) ! Total (global) number of PW in each direction
integer, intent(in) :: nsub(:) ! Number of subdomains in each direction
! Temporary input array, will get trashed by pfft3
complex(dp) :: tmp(size(x,1), size(x,2), size(x,3))
tmp = xG
! Calculates sum_{k=0}^{N-1} e^{2*pi*i*k*n/N}*X(k)
call pifft3(tmp, x, commy, commz, Ng, nsub)
! The result is already normalized
end subroutine

end module
