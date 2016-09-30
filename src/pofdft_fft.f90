module pofdft_fft
use types, only: dp
use constants, only: pi, Ha2eV
use ffte, only: dp_ffte
use pffte, only: pfft3_init, pfft3, pifft3
use utils, only: assert, stop_error, clock, allocate_mold
use xc, only: xc_pz
use ofdft, only: f
use ofdft_fft, only: update_fletcher_reeves, update_polak_ribiere
use optimize, only: bracket, brent, parabola_vertex
use mpi2, only: MPI_DOUBLE_PRECISION, mpi_allreduce, MPI_SUM, MPI_MAX
use integration, only: integrate_trapz_1
implicit none
private
public pfft3_init, preal2fourier, pfourier2real, real_space_vectors, &
    reciprocal_space_vectors, calculate_myxyz, pintegral, pintegralG, &
    free_energy, free_energy_min, radial_potential_fourier, psum, pmaxval

interface preal2fourier
    module procedure preal2fourier_real
    module procedure preal2fourier_complex
end interface

interface pfourier2real
    module procedure pfourier2real_real
    module procedure pfourier2real_complex
end interface

integer :: fft_counter
real(dp) :: fft_time = 0

contains

subroutine preal2fourier_complex(x, xG, commy, commz, Ng, nsub)
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
real(dp) :: t1, t2
t1 = clock()
tmp = x
! Calculates sum_{n=0}^{N-1} e^{-2*pi*i*k*n/N}*x(n)
call pfft3(tmp, xG, commy, commz, Ng, nsub)
xG = xG / product(Ng)     ! The proper normalization is to divide by N
t2 = clock()
fft_counter = fft_counter + 1
fft_time = fft_time + t2-t1
end subroutine

subroutine preal2fourier_real(x, xG, commy, commz, Ng, nsub)
! The same as preal2fourier_complex, but discard the imaginary part of "x"
real(dp), intent(in) :: x(:, :, :)
complex(dp), intent(out) :: xG(:, :, :)
integer, intent(in) :: commy, commz ! communicators in y, z directions
integer, intent(in) :: Ng(:) ! Total (global) number of PW in each direction
integer, intent(in) :: nsub(:) ! Number of subdomains in each direction
! Temporary input array, will get trashed by pfft3
complex(dp) :: tmp(size(x,1), size(x,2), size(x,3))
real(dp) :: t1, t2
t1 = clock()
tmp = x
! Calculates sum_{n=0}^{N-1} e^{-2*pi*i*k*n/N}*x(n)
call pfft3(tmp, xG, commy, commz, Ng, nsub)
xG = xG / product(Ng)     ! The proper normalization is to divide by N
t2 = clock()
fft_counter = fft_counter + 1
fft_time = fft_time + t2-t1
end subroutine

subroutine pfourier2real_complex(xG, x, commy, commz, Ng, nsub)
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
real(dp) :: t1, t2
t1 = clock()
tmp = xG
! Calculates sum_{k=0}^{N-1} e^{2*pi*i*k*n/N}*X(k)
call pifft3(tmp, x, commy, commz, Ng, nsub)
! The result is already normalized
t2 = clock()
fft_counter = fft_counter + 1
fft_time = fft_time + t2-t1
end subroutine

subroutine pfourier2real_real(xG, x, commy, commz, Ng, nsub)
! The same as pfourier2real_complex, but discard the imaginary part of "x"
complex(dp), intent(in) :: xG(:, :, :)
real(dp), intent(out) :: x(:, :, :)
integer, intent(in) :: commy, commz ! communicators in y, z directions
integer, intent(in) :: Ng(:) ! Total (global) number of PW in each direction
integer, intent(in) :: nsub(:) ! Number of subdomains in each direction
! Temporary input array, will get trashed by pfft3
complex(dp) :: tmp(size(x,1), size(x,2), size(x,3))
call pfourier2real_complex(xG, tmp, commy, commz, Ng, nsub)
x = real(tmp, dp)
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

real(dp) function psum(comm, f) result(r)
! Calculates the sum over 'f' in parallel, returns the answer on all
! processors.
integer, intent(in) :: comm
real(dp), intent(in) :: f(:,:,:)
real(dp) :: myr
integer :: ierr
myr = sum(f)
call mpi_allreduce(myr, r, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
end function

real(dp) function pmaxval(comm, f) result(r)
! Calculates the maxval over 'f' in parallel, returns the answer on all
! processors.
integer, intent(in) :: comm
real(dp), intent(in) :: f(:,:,:)
real(dp) :: myr
integer :: ierr
myr = maxval(f)
call mpi_allreduce(myr, r, 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm, ierr)
end function

real(dp) function pintegral(comm, L, f, Ng) result(r)
! Calculates the integral over 'f' in parallel, returns the answer on all
! processors.
integer, intent(in) :: comm
real(dp), intent(in) :: L(:), f(:,:,:)
integer, intent(in) :: Ng(:)
r = psum(comm, f) * product(L/Ng)
end function

real(dp) function pintegralG(comm, L, fG) result(r)
! Calculates the integral over 'fG' in reciprocal space in parallel, returns
! the answer on all processors.
integer, intent(in) :: comm
real(dp), intent(in) :: L(:)
real(dp), intent(in) :: fG(:, :, :)
r = psum(comm, fG) * product(L)
end function

real(dp) pure function integrate(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)
s = integrate_trapz_1(Rp, f)
end function

subroutine radial_potential_fourier(R, V, L, Z, Ng, myxyz, VenG, V0)
! Takes radial potential (given at origin) defined by:
!   V(R) on a grid R
!   Z/r for r > maxval(R)
! and calculates a 3D Fourier transform of it into the VenG grid. 'L' is the
! length of the box.
real(dp), intent(in) :: R(:), V(:), L(:), Z
integer, intent(in) :: Ng(:), myxyz(:)
real(dp), intent(out) :: VenG(:, :, :)
real(dp), intent(out) :: V0
integer :: Ng_local(3), ijk_global(3), i, j, k, idx
real(dp) :: Rp(size(R)), dk(3), Rc, w, Vk(0:3*(maxval(Ng)/2+1)**2)
! Rp is the derivative of the mesh R'(t), which for uniform mesh is equal to
! the mesh step (rmax-rmin)/N:
Rp = (R(size(R)) - R(1)) / (size(R)-1)
Rc = R(size(R))
Ng_local = shape(VenG)
dk = 2*pi/L

! V0 =      \int V(r) - Z/r d^3x
!    = 4*pi*\int_0^Rc (V(r) - Z/r) r^2 dr
!    = 4*pi*\int_0^Rc V(r) r^2 dr - Z*Rc^2 / 2
! So these two are equivalent, we use the one that is faster:
!V0 = -4*pi*integrate(Rp, R**2*(V-Z/R))
V0 = -4*pi*(integrate(Rp, R**2*V) - Z*Rc**2/2)

! We prepare the values of the radial Fourier transform on a 3D grid
Vk(0) = 0
do i = 1, 3*(maxval(Ng)/2+1)**2
    ! TODO: this formula recovers the sqrt(i)*dk result when dk is a constant
    ! vector. Figure out how to extend this algorithm when dk (or L) is not a
    ! constant vector.
    w = sqrt(i*sum(dk**2)/3)
    Vk(i) = 4*pi/(product(L) * w**2) * (w*integrate(Rp, R*sin(w*R)*V) + Z*cos(w*Rc))
end do

! We fill out the 3D grid using the values from Vk
do k = 1, Ng_local(3)
do j = 1, Ng_local(2)
do i = 1, Ng_local(1)
    ijk_global = [i, j, k] + myxyz*Ng_local
    idx = sum((ijk_global - 1 - Ng*nint((ijk_global-1.5_dp)/Ng))**2)
    VenG(i, j, k) = Vk(idx)
end do
end do
end do
end subroutine

subroutine poisson_kernel(myid, n, neG, G2, VeeG)
! Calculates VeeG = 4*pi*neG / G2, but skips the (1,1,1) term on zero rank
integer, intent(in) :: myid, n
complex(dp), intent(in) :: neG(n)
real(dp), intent(in) :: G2(n)
complex(dp), intent(out) :: VeeG(n)
if (myid == 0) then
    VeeG(1) = 0
    VeeG(2:) = 4*pi*neG(2:) / G2(2:)
else
    VeeG = 4*pi*neG / G2
end if
end subroutine

subroutine free_energy(myid, comm, commy, commz, Ng, nsub, &
        L, G2, T_au, VenG, ne, Eee, Een, Ts, Exc, Etot, dFdn, &
        calc_value, calc_derivative, vW, lambda1, lambda2, EvW)
use ofdft, only: f
integer, intent(in) :: myid, comm, commy, commz, Ng(:), nsub(:)
real(dp), intent(in) :: L(:), G2(:, :, :), T_au, ne(:, :, :)
complex(dp), intent(in) :: VenG(:, :, :)
real(dp), intent(out) :: Eee, Een, Ts, Exc, Etot
logical, intent(in) :: calc_value, calc_derivative
! include von Weizsäcker term, default .false.
logical, intent(in), optional :: vW
! The coefficient of the von Weizsäcker term. lambda1 is used when calculating
! energy, lambda2 when calculating the potential
real(dp), intent(in), optional :: lambda1, lambda2
! If vW == .true., this contains the von Weizsäcker part of the energy
real(dp), intent(out), optional :: EvW

! dFdn returns "delta F / delta n", the functional derivative with respect to
! the density "n". Use the relation
!     d/dpsi = 2 psi d/dn
! to obtain the derivative with respect to psi (i.e. multiply Hpsi by 2*psi).
real(dp), intent(out) :: dFdn(:, :, :)

real(dp), dimension(size(VenG,1), size(VenG,2), size(VenG,3)) :: y, F0, &
    exc_density, Ven_ee, Vxc, dF0dn, psi, d2psi, dFvWdn
complex(dp), dimension(size(VenG,1), size(VenG,2), size(VenG,3)) :: &
    neG, VeeG, psiG
real(dp) :: beta, dydn
logical :: vW_
vW_ = .false.
if (present(vW)) vW_ = vW
if (vW_) then
    call assert(present(lambda1))
    call assert(present(lambda2))
    call assert(present(EvW))
end if

call assert(calc_value .or. calc_derivative)

call preal2fourier(ne, neG, commy, commz, Ng, nsub)

call poisson_kernel(myid, size(neG), neG, G2, VeeG)

if (vW_) then
    psi = sqrt(ne)
    call preal2fourier(psi, psiG, commy, commz, Ng, nsub)
end if

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
    if (vW_) then
        EvW = pintegralG(comm, L, G2*abs(psiG)**2)/2 * lambda1
        Etot = Etot + EvW
    end if
end if

if (calc_derivative) then
    ! Calculate the derivative
    dydn = pi**2 / sqrt(2._dp) * beta**(3._dp/2)
    ! F0 = ne / beta * f(y)
    ! d F0 / d n =
    dF0dn = 1 / beta * f(y) + ne / beta * f(y, deriv=.true.) * dydn

    call pfourier2real(VenG+VeeG, Ven_ee, commy, commz, Ng, nsub)

    dFdn = dF0dn + Ven_ee + Vxc

    if (vW_) then
        call pfourier2real(-G2*psiG, d2psi, commy, commz, Ng, nsub)
        dFvWdn = -d2psi/(2*psi) * lambda2
        dFdn = dFdn + dFvWdn
    end if
end if
end subroutine

subroutine precalc_C1C2C2(myid, G2, psi, eta, C1, C2, C3, commy, commz, Ng, &
        nsub)
integer, intent(in) :: myid
real(dp), dimension(:, :, :), intent(in) :: G2, psi, eta
real(dp), dimension(:, :, :), intent(out) :: C1, C2, C3
integer, intent(in) :: commy, commz, Ng(:), nsub(:)
complex(dp), dimension(size(psi,1), size(psi,2), size(psi,3)) :: neG, VeeG
call preal2fourier(psi**2, neG, commy, commz, Ng, nsub)
call poisson_kernel(myid, size(neG), neG, G2, VeeG)
call pfourier2real(VeeG, C1, commy, commz, Ng, nsub)

call preal2fourier(psi*eta, neG, commy, commz, Ng, nsub)
call poisson_kernel(myid, size(neG), neG, G2, VeeG)
call pfourier2real(VeeG, C2, commy, commz, Ng, nsub)

call preal2fourier(eta**2, neG, commy, commz, Ng, nsub)
call poisson_kernel(myid, size(neG), neG, G2, VeeG)
call pfourier2real(VeeG, C3, commy, commz, Ng, nsub)
end subroutine

subroutine calc_F(comm, Ng, L, T_au, Ven, C1, C2, C3, psi, eta, theta, F_)
! Calculates F[psi(theta)] efficiently (no FFT needed, just integrals)
! See the function find_theta() for an example how to prepare the C1, C2 and C3
! constants and how to use this subroutine.
integer, intent(in) :: comm, Ng(:)
real(dp), intent(in) :: L(:), T_au, psi(:, :, :), eta(:, :, :)
real(dp), intent(in) :: C1(:, :, :), C2(:, :, :), C3(:, :, :)
real(dp), intent(in) :: theta
real(dp), intent(in) :: Ven(:, :, :)
real(dp), intent(out) :: F_
real(dp), dimension(size(psi, 1), size(psi, 2), size(psi, 3)) :: psi_, ne, y, &
    Vee, exc_density, Vxc
real(dp) :: beta
psi_ = psi*cos(theta)+eta*sin(theta)
ne = psi_**2

call xc_pz(ne, exc_density, Vxc)
beta = 1/T_au
y = pi**2 / sqrt(2._dp) * beta**(3._dp/2) * ne
if (any(y < 0)) call stop_error("Density must be positive")
Vee = C1*cos(theta)**2 + C2*sin(2*theta) + C3*sin(theta)**2

F_ = pintegral(comm, L, ne / beta * f(y) + (Ven+Vee/2+exc_density)*ne, Ng)
end subroutine

subroutine free_energy_min(myid, comm, commy, commz, Ng, nsub, &
        Nelec, Natom, L, G2, T_au, VenG, ne, energy_eps, &
        Eee, Een, Ts, Exc, Etot, cg_iter, vW, lambda, EvW)
! Minimize the electronic free energy using the initial condition 'ne'. Returns
! the ground state in 'ne'. The free energy is returned in Etot, and it's
! components are returned in Eee, Een, Ts and Exc. The relation is:
!
!     Etot = Ts + Een + Eee + Exc
!
real(dp), intent(in) :: Nelec ! Number of electrons
integer, intent(in) :: Natom ! Number of atoms
real(dp), intent(in) :: L(:), G2(:, :, :), T_au, energy_eps
real(dp), intent(inout) :: ne(:, :, :)
complex(dp), intent(in) :: VenG(:, :, :)
real(dp), intent(out) :: Eee, Een, Ts, Exc, Etot
integer, intent(out) :: cg_iter ! # of CG iterations needed to converge
integer, intent(in) :: myid, comm, commy, commz, Ng(:), nsub(:)
! include von Weizsäcker term, default .false.
logical, intent(in), optional :: vW
! The coefficient of the von Weizsäcker term
real(dp), intent(in), optional :: lambda
! If vW == .true., this contains the von Weizsäcker part of the energy
real(dp), intent(out), optional :: EvW

real(dp), allocatable :: free_energies(:)
real(dp), allocatable, dimension(:, :, :) :: Hpsi, &
    psi, psi_, psi_prev, ksi, ksi_prev, phi, phi_prime, eta
integer :: iter, max_iter
real(dp) :: mu, last2, last3, brent_eps, free_energy_, &
    gamma_d, gamma_n, theta, theta_a, theta_b, theta_c, fa, fb, fc
real(dp) :: f2
real(dp) :: psi_norm
integer :: update_type
!real(dp) :: A, B
real(dp) :: t1, t2, t11, t12
real(dp), dimension(size(VenG,1), size(VenG,2), size(VenG,3)) :: Ven, C1, C2, C3
logical :: func_use_fft
logical :: vW_
vW_ = .false.
if (present(vW)) vW_ = vW

brent_eps = 1e-3_dp
max_iter = 2000
update_type = update_polak_ribiere
func_use_fft = .true.

last2 = 0
last3 = 0

call allocate_mold(Hpsi, ne)
call allocate_mold(psi, ne)
call allocate_mold(psi_, ne)
call allocate_mold(psi_prev, ne)
call allocate_mold(phi, ne)
call allocate_mold(phi_prime, ne)
call allocate_mold(ksi, ne)
call allocate_mold(ksi_prev, ne)
call allocate_mold(eta, ne)

psi = sqrt(ne)
psi_norm = pintegral(comm, L, psi**2, Ng)
if (myid == 0) print *, "Initial norm of psi:", psi_norm
psi = sqrt(Nelec / psi_norm) * psi
psi_norm = pintegral(comm, L, psi**2, Ng)
if (myid == 0) print *, "norm of psi:", psi_norm
! This returns H[n] = delta F / delta n, we save it to the Hpsi variable to
! save space:
call free_energy(myid, comm, commy, commz, Ng, nsub, &
    L, G2, T_au, VenG, psi**2, Eee, Een, Ts, Exc, free_energy_, &
    Hpsi, calc_value=.false., calc_derivative=.true., vW=vW_, &
    lambda1=lambda, lambda2=lambda, EvW=EvW)
! Hpsi = H[psi] = delta F / delta psi = 2*H[n]*psi, due to d/dpsi = 2 psi d/dn
Hpsi = Hpsi * 2*psi
mu = 1._dp / Nelec * pintegral(comm, L, 0.5_dp * psi * Hpsi, Ng)
ksi = 2*mu*psi - Hpsi
phi = ksi
phi_prime = phi - 1._dp / Nelec *  pintegral(comm, L, phi * psi, Ng) * psi
eta = sqrt(Nelec / pintegral(comm, L, phi_prime**2, Ng)) * phi_prime
theta = pi/2
!print *, "Summary of energies [a.u.]:"
!print "('    Ts   = ', f14.8)", Ts
!print "('    Een  = ', f14.8)", Een
!print "('    Eee  = ', f14.8)", Eee
!print "('    Exc  = ', f14.8)", Exc
!print *, "   ---------------------"
!print "('    Etot = ', f14.8, ' a.u.')", free_energy_
allocate(free_energies(max_iter))
gamma_n = 0
do iter = 1, max_iter
    t1 = clock()
    fft_counter = 0
    fft_time = 0
    !Hpsi = Hpsi/(2*psi)
    ! Formula (36) in Jiang & Yang
    !A = integral(L, psi*Hpsi*psi) - integral(L, eta*Hpsi*eta)
    !B = 2*integral(L, eta*Hpsi*psi)
    !print *, "theta? =", 0.5_dp * atan(B/A)
    theta_a = 0
    theta_b = mod(theta, 2*pi)
    t11 = clock()
    if (iter == 1) then
        if (func_use_fft) then
            call bracket(func_fft, theta_a, theta_b, theta_c, fa, fb, fc, 100._dp, 20, verbose=.false.)
            call brent(func_fft, theta_a, theta_b, theta_c, brent_eps, 50, theta, &
            free_energy_, verbose=.false.)
        else
            call pfourier2real(VenG, Ven, commy, commz, Ng, nsub)
            call precalc_C1C2C2(myid, G2, psi, eta, C1, C2, C3, &
                commy, commz, Ng, nsub)
            call bracket(func_nofft, theta_a, theta_b, theta_c, fa, fb, fc, 100._dp, 20, verbose=.false.)
            call brent(func_nofft, theta_a, theta_b, theta_c, brent_eps, 50, theta, &
            free_energy_, verbose=.false.)
        end if
    else
        call bracket(func_fft, theta_a, theta_b, theta_c, fa, fb, fc, 100._dp, 20, verbose=.false.)
        call parabola_vertex(theta_a, fa, theta_b, fb, theta_c, fc, theta, f2)
    end if
    t12 = clock()
    ! TODO: We probably don't need to recalculate free_energy_ here:
    psi_prev = psi
    psi = cos(theta) * psi + sin(theta) * eta
    call free_energy(myid, comm, commy, commz, Ng, nsub, &
        L, G2, T_au, VenG, psi**2, Eee, Een, Ts, Exc, &
        free_energy_, Hpsi, calc_value=.true., calc_derivative=.true., &
        vW=vW_, lambda1=lambda, lambda2=lambda, EvW=EvW)
!    print *, "Iteration:", iter
!    psi_norm = integral(L, psi**2)
!    print *, "Norm of psi:", psi_norm
!    print *, "mu =", mu
!    print *, "|ksi| =", sqrt(gamma_n)
!    print *, "theta =", theta
!    print *, "Summary of energies [a.u.]:"
!    print "('    Ts   = ', f14.8)", Ts
!    print "('    Een  = ', f14.8)", Een
!    print "('    Eee  = ', f14.8)", Eee
!    print "('    Exc  = ', f14.8)", Exc
!    print *, "   ---------------------"
    free_energies(iter) = free_energy_ / Natom
    if (iter > 1) then
        last2 = maxval(free_energies(iter-1:iter)) - &
            minval(free_energies(iter-1:iter))
    end if
    if (iter > 2) then
        last3 = maxval(free_energies(iter-2:iter)) - &
            minval(free_energies(iter-2:iter))
    end if
    if (myid == 0) then
        print "('# ', i3, ' Etot/atom = ', f18.8, ' eV; last2 = ', es10.2, ' last3 = ',es10.2)", &
        iter, free_energy_ * Ha2eV / Natom, last2 * Ha2eV, last3 * Ha2eV
    end if
    if (iter > 3) then
        if (last3 < energy_eps) then
            ne = psi**2
            Etot = free_energy_
            cg_iter = iter
            return
        end if
    end if
    Hpsi = Hpsi * 2*psi ! d/dpsi = 2 psi d/dn
    mu = 1._dp / Nelec * pintegral(comm, L, 0.5_dp * psi * Hpsi, Ng)
    ksi_prev = ksi
    ksi = 2*mu*psi - Hpsi
    select case(update_type)
        case(update_fletcher_reeves) ! Fletcher-Reeves
            gamma_n = pintegral(comm, L, ksi**2, Ng)
        case(update_polak_ribiere)   ! Polak-Ribiere
            gamma_n = max(pintegral(comm, L, ksi*(ksi-ksi_prev), Ng), 0._dp)
        case default
            call stop_error("Unknown update type.")
    end select
    gamma_d = pintegral(comm, L, ksi_prev**2, Ng)
    phi = ksi + gamma_n / gamma_d * phi
    phi_prime = phi - 1._dp / Nelec * pintegral(comm, L, phi * psi, Ng) * psi
    eta = sqrt(Nelec / pintegral(comm, L, phi_prime**2, Ng)) * phi_prime
    t2 = clock()
    if (myid == 0) then
        print "('time: ', f6.3, ' fft: ', f6.3, ' (', f5.2, '%) # fft: ', i3, ' line search', f6.3)", &
            t2-t1, fft_time, fft_time / (t2-t1) * 100, fft_counter, t12-t11
    end if
end do
call stop_error("free_energy_minimization: The maximum number of iterations exceeded.")

contains

    real(dp) function func_nofft(theta) result(energy)
    real(dp), intent(in) :: theta
    call calc_F(comm, Ng, L, T_au, Ven, C1, C2, C3, psi, eta, theta, energy)
    end function

    real(dp) function func_fft(theta) result(energy)
    real(dp), intent(in) :: theta
    psi_ = cos(theta) * psi + sin(theta) * eta
    call free_energy(myid, comm, commy, commz, Ng, nsub, &
        L, G2, T_au, VenG, psi_**2, Eee, Een, Ts, Exc, &
        energy, Hpsi, calc_value=.true., calc_derivative=.false., &
        vW=vW_, lambda1=lambda, lambda2=lambda, EvW=EvW)
    end function

end subroutine

end module
