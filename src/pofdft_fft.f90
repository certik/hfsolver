module pofdft_fft
use types, only: dp
use constants, only: pi
use ffte, only: dp_ffte
use pffte, only: pfft3_init, pfft3, pifft3
implicit none
private
public pfft3_init, preal2fourier, pfourier2real, real_space_vectors, &
    reciprocal_space_vectors


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

subroutine real_space_vectors_3d(L, X)
! Calculates the real space vectors in the box [0, L]^3
real(dp), intent(in) :: L(:)
! X(i, j, k, :) is the 3D position vector of the point with index (i, j, k)
real(dp), intent(out) :: X(:, :, :, :)
integer :: Ng(3), i, j, k
Ng = [size(X,1), size(X,2), size(X,3)]
forall(i=1:Ng(1), j=1:Ng(2), k=1:Ng(3))
    X(i, j, k, :) = ([i, j, k] - 1) * L / Ng
end forall
end subroutine

subroutine reciprocal_space_vectors_3d(L, G, G2)
real(dp), intent(in) :: L(:)
! G(:, :, :, i) where i=1, 2, 3 are the x, y, z components
! G2(:, :, :) are the squares of G
real(dp), intent(out) :: G(:, :, :, :), G2(:, :, :)
integer :: Ng(3), i, j, k
Ng = shape(G2)
forall(i=1:Ng(1), j=1:Ng(2), k=1:Ng(3))
    G(i, j, k, :) = 2*pi/L * ([i,j,k] - 1 - Ng*nint(([i,j,k]-1.5_dp)/Ng))
end forall
G(1, 1, 1, :) = 1 ! To avoid division by 0
G2 = G(:,:,:,1)**2 + G(:,:,:,2)**2 + G(:,:,:,3)**2
end subroutine

subroutine real_space_vectors(L, X, Ng, myxyz)
real(dp), intent(in) :: L(:)
real(dp), intent(out) :: X(:, :, :, :)
integer, intent(in) :: Ng(:), myxyz(:)
integer:: Ng_local(3), ijk_global(3), i, j, k
real(dp) :: X_global(Ng(1), Ng(2), Ng(3), 3)
call real_space_vectors_3d(L, X_global)
Ng_local = [size(X,1), size(X,2), size(X,3)]
do k = 1, size(X, 3)
do j = 1, size(X, 2)
do i = 1, size(X, 1)
    ijk_global = [i, j, k] + myxyz*Ng_local
    X(i,j,k,:) = X_global(ijk_global(1), ijk_global(2), ijk_global(3), :)
end do
end do
end do
end subroutine

subroutine reciprocal_space_vectors(L, G, G2, Ng, myxyz)
real(dp), intent(in) :: L(:)
! G(:, :, :, i) where i=1, 2, 3 are the x, y, z components
! G2(:, :, :) are the squares of G
real(dp), intent(out) :: G(:, :, :, :), G2(:, :, :)
integer, intent(in) :: Ng(:), myxyz(:)
integer:: Ng_local(3), ijk_global(3), i, j, k
real(dp) :: G_global(Ng(1), Ng(2), Ng(3), 3), G2_global(Ng(1),Ng(2),Ng(3))
call reciprocal_space_vectors_3d(L, G_global, G2_global)
Ng_local = [size(G,1), size(G,2), size(G,3)]
do k = 1, size(G, 3)
do j = 1, size(G, 2)
do i = 1, size(G, 1)
    ijk_global = [i, j, k] + myxyz*Ng_local
    G(i,j,k,:) = G_global(ijk_global(1), ijk_global(2), ijk_global(3), :)
    G2(i,j,k) = G2_global(ijk_global(1), ijk_global(2), ijk_global(3))
end do
end do
end do
end subroutine

end module
