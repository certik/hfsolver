module pofdft_fft
use types, only: dp
use constants, only: pi
use ffte, only: dp_ffte
use pffte, only: pfft3_init, pfft3, pifft3
use utils, only: assert, stop_error
use xc, only: xc_pz
use mpi2, only: MPI_DOUBLE_PRECISION, mpi_allreduce, MPI_SUM
implicit none
private
public pfft3_init, preal2fourier, pfourier2real, real_space_vectors, &
    reciprocal_space_vectors, calculate_myxyz, pintegral, pintegralG


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
! normalization as defined by preal2fourier(), i.e. in the following calls, 'x2'
! will be equal to 'x':
! call preal2fourier(x, xG)
! call pfourier2real(xG, x2)
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

function calculate_myxyz(myid, nsub) result(myxyz)
integer, intent(in) :: myid, nsub(3)
integer :: myxyz(3)
myxyz = [0, mod(myid, nsub(2)), myid/nsub(2)]
end function

subroutine real_space_vectors(L, X, Ng, myxyz)
real(dp), intent(in) :: L(:)
real(dp), intent(out) :: X(:, :, :, :)
integer, intent(in) :: Ng(:), myxyz(:)
integer:: Ng_local(3), ijk_global(3), i, j, k
Ng_local = [size(X,1), size(X,2), size(X,3)]
do k = 1, size(X, 3)
do j = 1, size(X, 2)
do i = 1, size(X, 1)
    ijk_global = [i, j, k] + myxyz*Ng_local
    X(i,j,k,:) = (ijk_global-1) * L / Ng
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
Ng_local = [size(G,1), size(G,2), size(G,3)]
do k = 1, size(G, 3)
do j = 1, size(G, 2)
do i = 1, size(G, 1)
    ijk_global = [i, j, k] + myxyz*Ng_local
    G(i,j,k,:) = 2*pi/L * (ijk_global - 1 - Ng*nint((ijk_global-1.5_dp)/Ng))
end do
end do
end do
G2 = G(:,:,:,1)**2 + G(:,:,:,2)**2 + G(:,:,:,3)**2
end subroutine

real(dp) function pintegral(comm, L, f, Ng) result(r)
! Calculates the integral over 'f' in parallel, returns the answer on all
! processors.
integer, intent(in) :: comm
real(dp), intent(in) :: L(:), f(:,:,:)
integer, intent(in) :: Ng(:)
real(dp) :: myr
integer :: ierr
myr = sum(f)
call mpi_allreduce(myr, r, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
r = r * product(L/Ng)
end function

real(dp) function pintegralG(comm, L, fG) result(r)
! Calculates the integral over 'fG' in reciprocal space in parallel, returns
! the answer on all processors.
integer, intent(in) :: comm
real(dp), intent(in) :: L(:)
real(dp), intent(in) :: fG(:, :, :)
real(dp) :: myr
integer :: ierr
!myr = sum(real(fG, dp))
myr = sum(fG)
call mpi_allreduce(myr, r, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
r = r * product(L)
end function

subroutine free_energy(comm, commy, commz, Ng, nsub, &
        L, G2, T_au, VenG, ne, Eee, Een, Ts, Exc, Etot, dFdn, &
        calc_value, calc_derivative)
use ofdft, only: f
integer, intent(in) :: comm, commy, commz, Ng(:), nsub(:)
real(dp), intent(in) :: L(:), G2(:, :, :), T_au, ne(:, :, :)
complex(dp), intent(in) :: VenG(:, :, :)
real(dp), intent(out) :: Eee, Een, Ts, Exc, Etot
logical, intent(in) :: calc_value, calc_derivative

! dFdn returns "delta F / delta n", the functional derivative with respect to
! the density "n". Use the relation
!     d/dpsi = 2 psi d/dn
! to obtain the derivative with respect to psi (i.e. multiply Hpsi by 2*psi).
real(dp), intent(out) :: dFdn(:, :, :)

real(dp), dimension(size(VenG,1), size(VenG,2), size(VenG,3)) :: y, F0, &
    exc_density, Ven_ee, Vxc, dF0dn
complex(dp), dimension(size(VenG,1), size(VenG,2), size(VenG,3)) :: &
    neG, VeeG, ne_, Ven_ee_
real(dp) :: beta, dydn
call assert(calc_value .or. calc_derivative)

ne_ = ne
call preal2fourier(ne_, neG, commy, commz, Ng, nsub)

VeeG = 4*pi*neG / G2
VeeG(1, 1, 1) = 0

beta = 1/T_au
! The density must be positive, the f(y) fails for negative "y". Thus we use
! ne.
y = pi**2 / sqrt(2._dp) * beta**(3._dp/2) * ne
if (any(y < 0)) call stop_error("Density must be positive")

call xc_pz(ne, exc_density, Vxc)

if (calc_value) then
    ! Hartree energy
    !Eee = integralG2(L, real(VeeG)*real(neG)+aimag(VeeG)*aimag(neG)) / 2
    Eee = pintegralG(comm, L, real(VeeG*conjg(neG), dp)) / 2
    ! Electron-nucleus energy
    !Een = integralG2(L, real(VenG)*real(neG)+aimag(VenG)*aimag(neG))
    Een = pintegralG(comm, L, real(VenG*conjg(neG), dp))

    ! Kinetic energy using Perrot parametrization
    F0 = ne / beta * f(y)
    Ts = pintegral(comm, L, F0, Ng)
    ! Exchange and correlation potential
    Exc = pintegral(comm, L, exc_density * ne, Ng)
    Etot = Ts + Een + Eee + Exc
end if

if (calc_derivative) then
    ! Calculate the derivative
    dydn = pi**2 / sqrt(2._dp) * beta**(3._dp/2)
    ! F0 = ne / beta * f(y)
    ! d F0 / d n =
    dF0dn = 1 / beta * f(y) + ne / beta * f(y, deriv=.true.) * dydn

    call pfourier2real(VenG+VeeG, Ven_ee_, commy, commz, Ng, nsub)
    Ven_ee = real(Ven_ee_, dp)

    dFdn = dF0dn + Ven_ee + Vxc
    dFdn = dFdn ! Do not multiply by L**3
end if
end subroutine

end module
