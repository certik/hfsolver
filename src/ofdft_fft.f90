module ofdft_fft

! Routines for orbital free density functional theory, discretization is
! done using FFT.
!
! We use 3D reciprocal grid (:, :, :). Variables in reciprocal (Fourier) space
! have G appended to them. E.g. Ven is the electron-nucleus potential in real
! space, VenG is the same potential in reciprocal space.

use types, only: dp
use constants, only: pi, i_, Ha2eV
use optimize, only: bracket, brent, parabola_vertex
use ofdft, only: f
use xc, only: xc_pz, xc_pz2
use utils, only: stop_error, assert, clock
use fourier, only: fft3_inplace, ifft3_inplace, &
    fft1_inplace => fft_pass_inplace, ifft1_inplace
use integration, only: integrate_trapz_1
implicit none
private
public reciprocal_space_vectors, free_energy, free_energy_min, &
    radial_potential_fourier, real2fourier, fourier2real, integralG, &
    logging_info, integral, real_space_vectors, radial_density_fourier, &
    vtk_save, update_fletcher_reeves, update_polak_ribiere

! Update types for nonlinear conjugate gradient method:
integer, parameter :: update_fletcher_reeves = 1
integer, parameter :: update_polak_ribiere   = 2

logical :: logging_info = .true.

interface real_space_vectors
    module procedure real_space_vectors_1d
    module procedure real_space_vectors_3d
end interface

interface reciprocal_space_vectors
    module procedure reciprocal_space_vectors_1d
    module procedure reciprocal_space_vectors_3d
end interface

interface integral
    module procedure integral3d, integral1d
end interface

interface integralG
    module procedure integralG_complex, integralG_real
    module procedure integralG_real_1d
end interface

interface real2fourier
    module procedure real2fourier_complex, real2fourier_real
    module procedure real2fourier1d_complex, real2fourier1d_real
end interface

interface fourier2real
    module procedure fourier2real_complex, fourier2real_real
    module procedure fourier2real1d_complex, fourier2real1d_real
end interface

integer :: fft_counter
real(dp) :: fft_time

contains

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

subroutine real_space_vectors_1d(L, X)
! Calculates the real space vectors in the box [0, L]
real(dp), intent(in) :: L
! X(i) is the 1D position vector of the point with index (i)
real(dp), intent(out) :: X(:)
integer :: Ng, i
Ng = size(X, 1)
forall(i=1:Ng) X(i) = (i-1) * L / Ng
end subroutine

subroutine reciprocal_space_vectors_3d(L, G, G2)
real(dp), intent(in) :: L(:)
! G(:, :, :, i) where i=1, 2, 3 are the x, y, z components
! G2(:, :, :) are the squares of G
real(dp), intent(out) :: G(:, :, :, :), G2(:, :, :)
integer :: Ng(3), i, j, k
Ng = shape(G2)
! The FFT frequencies are ordered as follows
!
!     G1D(1:Ng) = 2*pi*k/L,
!
! where
!
!     k = 0, 1, ..., Ng/2-1, -Ng/2, -Ng/2+1, ..., -1 if Ng is even,
!     k = 0, 1, ..., (Ng-1)/2, -(Ng-1)/2, -(Ng-1)/2+1, ..., -1 if Ng is odd.
!
! E.g. for
!
!     Ng=8 we get k = 0, 1, 2, 3, -4, -3, -2, -1.
!     Ng=9 we get k = 0, 1, 2, 3, 4, -4, -3, -2, -1.
!
! 'k' can be generated using
!
!     forall(i=1:Ng) k(i) = i-1 - Ng*nint((i-0.5_dp)/Ng)    ! if Ng is even
!     forall(i=1:Ng) k(i) = i-1 - Ng*nint((i-1.5_dp)/Ng)    ! if Ng is odd
!
! Note that the Ng/2 frequency can go both ways: +/- Ng/2. Below we take it with
! the "+" sign:
!
!     k = 0, 1, ..., Ng/2-1, +Ng/2, -Ng/2+1, ..., -1 if n is even
!
! E.g. for
!
!     Ng=8 we get k = 0, 1, 2, 3, 4, -3, -2, -1.
!     Ng=9 we get k = 0, 1, 2, 3, 4, -4, -3, -2, -1.
!
! These can be generated using (this holds for 'Ng' both even and odd):
!
!     forall(i=1:Ng) k(i) = i-1 - Ng*nint((i-1.5_dp)/Ng)
!
! That is the formula that we use below:
forall(i=1:Ng(1), j=1:Ng(2), k=1:Ng(3))
    G(i, j, k, :) = 2*pi/L * ([i,j,k] - 1 - Ng*nint(([i,j,k]-1.5_dp)/Ng))
end forall
G2 = G(:,:,:,1)**2 + G(:,:,:,2)**2 + G(:,:,:,3)**2
end subroutine

subroutine reciprocal_space_vectors_1d(L, G, G2)
real(dp), intent(in) :: L
! G(:) are the x components
! G2(:) are the squares of G
real(dp), intent(out) :: G(:), G2(:)
integer :: Ng, i
Ng = size(G, 1)
forall(i=1:Ng) G(i) = 2*pi/L * (i-1-Ng*nint((i-1.5_dp)/Ng))
G2 = G**2
end subroutine

subroutine real2fourier_complex(x, xG)
! Calculates Discrete Fourier Transform, with the same normalization as the
! Fourier Transform for periodic systems (which is Fourier Series).
complex(dp), intent(in) :: x(:, :, :)
complex(dp), intent(out) :: xG(:, :, :)
integer :: Ng
real(dp) :: t1, t2
t1 = clock()
Ng = size(x, 1)
xG = x
call fft3_inplace(xG) ! Calculates sum_{n=0}^{N-1} e^{-2*pi*i*k*n/N}*x(n)
xG = xG / size(x)     ! The proper normalization is to divide by N
t2 = clock()
fft_counter = fft_counter + 1
fft_time = fft_time + t2-t1
end subroutine

subroutine real2fourier1d_complex(x, xG)
! Calculates Discrete Fourier Transform, with the same normalization as the
! Fourier Transform for periodic systems (which is Fourier Series).
complex(dp), intent(in) :: x(:)
complex(dp), intent(out) :: xG(:)
integer :: Ng
real(dp) :: t1, t2
t1 = clock()
Ng = size(x, 1)
xG = x
call fft1_inplace(xG) ! Calculates sum_{n=0}^{N-1} e^{-2*pi*i*k*n/N}*x(n)
xG = xG / size(x)     ! The proper normalization is to divide by N
t2 = clock()
fft_counter = fft_counter + 1
fft_time = fft_time + t2-t1
end subroutine

subroutine real2fourier_real(x, xG)
! Calculates Discrete Fourier Transform, with the same normalization as the
! Fourier Transform for periodic systems (which is Fourier Series).
real(dp), intent(in) :: x(:, :, :)
complex(dp), intent(out) :: xG(:, :, :)
integer :: Ng
real(dp) :: t1, t2
t1 = clock()
Ng = size(x, 1)
xG = x
call fft3_inplace(xG) ! Calculates sum_{n=0}^{N-1} e^{-2*pi*i*k*n/N}*x(n)
xG = xG / size(x)     ! The proper normalization is to divide by N
t2 = clock()
fft_counter = fft_counter + 1
fft_time = fft_time + t2-t1
end subroutine

subroutine real2fourier1d_real(x, xG)
! Calculates Discrete Fourier Transform, with the same normalization as the
! Fourier Transform for periodic systems (which is Fourier Series).
real(dp), intent(in) :: x(:)
complex(dp), intent(out) :: xG(:)
integer :: Ng
real(dp) :: t1, t2
t1 = clock()
Ng = size(x, 1)
xG = x
call fft1_inplace(xG) ! Calculates sum_{n=0}^{N-1} e^{-2*pi*i*k*n/N}*x(n)
xG = xG / size(x)     ! The proper normalization is to divide by N
t2 = clock()
fft_counter = fft_counter + 1
fft_time = fft_time + t2-t1
end subroutine

subroutine fourier2real_complex(xG, x)
! Calculates Inverse Discrete Fourier Transform. xG must follow the same
! normalization as defined by real2fourier(), i.e. in the following calls, 'x2'
! will be equal to 'x':
! call real2fourier(x, xG)
! call fourier2real(xG, x2)
complex(dp), intent(in) :: xG(:, :, :)
complex(dp), intent(out) :: x(:, :, :)
real(dp) :: t1, t2
t1 = clock()
x = xG
call ifft3_inplace(x) ! Calculates sum_{k=0}^{N-1} e^{2*pi*i*k*n/N}*X(k)
! The result is already normalized
t2 = clock()
fft_counter = fft_counter + 1
fft_time = fft_time + t2-t1
end subroutine

subroutine fourier2real1d_complex(xG, x)
! Calculates Inverse Discrete Fourier Transform. xG must follow the same
! normalization as defined by real2fourier(), i.e. in the following calls, 'x2'
! will be equal to 'x':
! call real2fourier(x, xG)
! call fourier2real(xG, x2)
complex(dp), intent(in) :: xG(:)
complex(dp), intent(out) :: x(:)
real(dp) :: t1, t2
t1 = clock()
x = xG
call ifft1_inplace(x) ! Calculates sum_{k=0}^{N-1} e^{2*pi*i*k*n/N}*X(k)
! The result is already normalized
t2 = clock()
fft_counter = fft_counter + 1
fft_time = fft_time + t2-t1
end subroutine

subroutine fourier2real_real(xG, x)
! Calculates Inverse Discrete Fourier Transform. xG must follow the same
! normalization as defined by real2fourier(), i.e. in the following calls, 'x2'
! will be equal to 'x':
! call real2fourier(x, xG)
! call fourier2real(xG, x2)
complex(dp), intent(in) :: xG(:, :, :)
real(dp), intent(out) :: x(:, :, :)
complex(dp) :: tmp(size(xG,1), size(xG,2), size(xG,3))
real(dp) :: t1, t2
t1 = clock()
tmp = xG
call ifft3_inplace(tmp) ! Calculates sum_{k=0}^{N-1} e^{2*pi*i*k*n/N}*X(k)
x = real(tmp, dp)       ! The result is already normalized
if (maxval(abs(aimag(tmp))) > 1e-12_dp) then
    if (logging_info) then
        print *, "INFO: fourier2real() max imaginary part:", maxval(aimag(tmp))
    end if
end if
! TODO: make this strict:
!call assert(maxval(abs(aimag(tmp))) < 1e-1_dp)
t2 = clock()
fft_counter = fft_counter + 1
fft_time = fft_time + t2-t1
end subroutine

subroutine fourier2real1d_real(xG, x)
! Calculates Inverse Discrete Fourier Transform. xG must follow the same
! normalization as defined by real2fourier(), i.e. in the following calls, 'x2'
! will be equal to 'x':
! call real2fourier(x, xG)
! call fourier2real(xG, x2)
complex(dp), intent(in) :: xG(:)
real(dp), intent(out) :: x(:)
complex(dp) :: tmp(size(xG,1))
real(dp) :: t1, t2
t1 = clock()
tmp = xG
call ifft1_inplace(tmp) ! Calculates sum_{k=0}^{N-1} e^{2*pi*i*k*n/N}*X(k)
x = real(tmp, dp)       ! The result is already normalized
if (maxval(abs(aimag(tmp))) > 1e-12_dp) then
    if (logging_info) then
        print *, "INFO: fourier2real() max imaginary part:", maxval(aimag(tmp))
    end if
end if
! TODO: make this strict:
!call assert(maxval(abs(aimag(tmp))) < 1e-1_dp)
t2 = clock()
fft_counter = fft_counter + 1
fft_time = fft_time + t2-t1
end subroutine

real(dp) function integral3d(L, f) result(r)
real(dp), intent(in) :: L, f(:, :, :)
r = sum(f) * L**3 / size(f)
end function

real(dp) function integral1d(L, f) result(r)
real(dp), intent(in) :: L, f(:)
r = sum(f) * L / size(f)
end function

real(dp) function integralG_complex(L, fG) result(r)
real(dp), intent(in) :: L
complex(dp), intent(in) :: fG(:, :, :)
complex(dp) :: s
s = sum(fG) * L**3
r = real(s, dp)
if (abs(aimag(s)) > 1e-12_dp) then
    if (logging_info) then
        print *, "INFO: integralG() imaginary part:", aimag(s)
    end if
end if
!if (abs(aimag(s)) > 1e-5_dp) then
!    print *, "aimag(s) =", aimag(s)
!    call stop_error("integralG(): Complex part is not negligible.")
!end if
end function

real(dp) function integralG_real(fG, L) result(r)
real(dp), intent(in) :: L, fG(:, :, :)
r = sum(fG) * L**3
end function

real(dp) function integralG_real_1d(fG, L) result(r)
real(dp), intent(in) :: L, fG(:)
r = sum(fG) * L
end function

subroutine calc_dFdtheta(L, T_au, Ven, C1, C2, C3, psi, eta, theta, dFdtheta)
! Calculates dF[psi(theta)]/dtheta efficiently (no FFT needed, just integrals)
! See the function find_theta() for an example how to prepare the C1, C2 and C3
! constants and how to use this subroutine.
real(dp), intent(in) :: L, T_au, psi(:, :, :), eta(:, :, :)
real(dp), intent(in) :: C1(:, :, :), C2(:, :, :), C3(:, :, :)
real(dp), intent(in) :: theta
real(dp), intent(in) :: Ven(:, :, :)
real(dp), intent(out) :: dFdtheta
real(dp), dimension(size(psi, 1), size(psi, 2), size(psi, 3)) :: psi_, ne

real(dp), dimension(size(Ven,1), size(Ven,2), size(Ven,3)) :: y, &
    Vee, Vxc, dF0dn
real(dp) :: beta, dydn
psi_ = psi*cos(theta)+eta*sin(theta)
ne = psi_**2

call xc_pz2(ne, Vxc)
beta = 1/T_au
y = pi**2 / sqrt(2._dp) * beta**(3._dp/2) * ne
if (any(y < 0)) call stop_error("Density must be positive")
dydn = pi**2 / sqrt(2._dp) * beta**(3._dp/2)
dF0dn = 1 / beta * f(y) + ne / beta * f(y, deriv=.true.) * dydn
Vee = C1*cos(theta)**2 + C2*sin(2*theta) + C3*sin(theta)**2

dFdtheta = integral(L, (-psi*sin(theta)+eta*cos(theta)) * &
    (dF0dn + Vee + Ven + Vxc) * L**3 * psi_)
end subroutine

subroutine calc_F(L, T_au, Ven, C1, C2, C3, psi, eta, theta, F_)
! Calculates F[psi(theta)] efficiently (no FFT needed, just integrals)
! See the function find_theta() for an example how to prepare the C1, C2 and C3
! constants and how to use this subroutine.
real(dp), intent(in) :: L, T_au, psi(:, :, :), eta(:, :, :)
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

F_ = integral(L, ne / beta * f(y) + (Ven+Vee/2+exc_density)*ne)
end subroutine

subroutine precalc_C1C2C2(G2, psi, eta, C1, C2, C3)
real(dp), dimension(:, :, :), intent(in) :: G2, psi, eta
real(dp), dimension(:, :, :), intent(out) :: C1, C2, C3
complex(dp), dimension(size(psi,1), size(psi,2), size(psi,3)) :: neG, VeeG
call real2fourier(psi**2, neG)
VeeG = 4*pi*neG / G2
VeeG(1, 1, 1) = 0
call fourier2real(VeeG, C1)

call real2fourier(psi*eta, neG)
VeeG = 4*pi*neG / G2
VeeG(1, 1, 1) = 0
call fourier2real(VeeG, C2)

call real2fourier(eta**2, neG)
VeeG = 4*pi*neG / G2
VeeG(1, 1, 1) = 0
call fourier2real(VeeG, C3)
end subroutine

real(dp) function free_energy_theta(L, G2, T_au, VenG, psi, eta, theta) &
        result(F_)
! Shows how to use calc_F. Typically you want to precalculate C1, C2, C3 and
! store them.
real(dp), intent(in) :: L, G2(:, :, :), T_au, psi(:, :, :), eta(:, :, :)
complex(dp), intent(in) :: VenG(:, :, :)
real(dp), intent(in) :: theta
real(dp), dimension(size(VenG,1), size(VenG,2), size(VenG,3)) :: Ven, C1, C2, C3
call fourier2real(VenG, Ven)
call precalc_C1C2C2(G2, psi, eta, C1, C2, C3)
call calc_F(L, T_au, Ven, C1, C2, C3, psi, eta, theta, F_)
end function

subroutine find_theta(L, G2, T_au, VenG, psi, eta)
! It prints the values of dF/dtheta using efficient evaluation
real(dp), intent(in) :: L, G2(:, :, :), T_au, psi(:, :, :), eta(:, :, :)
complex(dp), intent(in) :: VenG(:, :, :)
real(dp) :: theta, dFdtheta
real(dp), dimension(size(VenG,1), size(VenG,2), size(VenG,3)) :: Ven, C1, C2, C3
integer :: i
real(dp) :: t1, t2
call fourier2real(VenG, Ven)
call precalc_C1C2C2(G2, psi, eta, C1, C2, C3)
print *, "dFdtheta:"
theta = 0
call cpu_time(t1)
do i = 1, 100
    call calc_dFdtheta(L, T_au, Ven, C1, C2, C3, psi, eta, theta, dFdtheta)
    print *, theta, dFdtheta
    theta = theta + pi/2/100
end do
call cpu_time(t2)
print *, "time: ", t2-t1
end subroutine

subroutine free_energy(L, G2, T_au, VenG, ne, Eee, Een, Ts, Exc, Etot, dFdn, &
    calc_value, calc_derivative)
real(dp), intent(in) :: L, G2(:, :, :), T_au, ne(:, :, :)
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
    neG, VeeG
real(dp) :: beta, dydn
call assert(calc_value .or. calc_derivative)

call real2fourier(ne, neG)

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
    Eee = integralG(L, VeeG*conjg(neG)) / 2
    ! Electron-nucleus energy
    !Een = integralG2(L, real(VenG)*real(neG)+aimag(VenG)*aimag(neG))
    Een = integralG(L, VenG*conjg(neG))

    ! Kinetic energy using Perrot parametrization
    F0 = ne / beta * f(y)
    Ts = integral(L, F0)
    ! Exchange and correlation potential
    Exc = integral(L, exc_density * ne)
    Etot = Ts + Een + Eee + Exc
end if

if (calc_derivative) then
    ! Calculate the derivative
    dydn = pi**2 / sqrt(2._dp) * beta**(3._dp/2)
    ! F0 = ne / beta * f(y)
    ! d F0 / d n =
    dF0dn = 1 / beta * f(y) + ne / beta * f(y, deriv=.true.) * dydn

    call fourier2real(VenG+VeeG, Ven_ee)

    dFdn = dF0dn + Ven_ee + Vxc
end if
end subroutine

subroutine radial_potential_fourier(R, V, L, Z, VenG, V0)
! Takes radial potential (given at origin) defined by:
!   V(R) on a grid R
!   Z/r for r > maxval(R)
! and calculates a 3D Fourier transform of it into the VenG grid. 'L' is the
! length of the box.
real(dp), intent(in) :: R(:), V(:), L, Z
real(dp), intent(out) :: VenG(:, :, :)
real(dp), intent(out) :: V0
integer :: Ng, i, j, k, idx
real(dp) :: Rp(size(R)), dk, Rc, w, Vk(0:3*(size(VenG, 1)/2+1)**2)
! Rp is the derivative of the mesh R'(t), which for uniform mesh is equal to
! the mesh step (rmax-rmin)/N:
Rp = (R(size(R)) - R(1)) / (size(R)-1)
Rc = R(size(R))
Ng = size(VenG, 1)
call assert(size(VenG, 2) == Ng)
call assert(size(VenG, 3) == Ng)
dk = 2*pi/L

! V0 =      \int V(r) - Z/r d^3x
!    = 4*pi*\int_0^Rc (V(r) - Z/r) r^2 dr
!    = 4*pi*\int_0^Rc V(r) r^2 dr - Z*Rc^2 / 2
! So these two are equivalent, we use the one that is faster:
!V0 = -4*pi*integrate(Rp, R**2*(V-Z/R))
V0 = -4*pi*(integrate(Rp, R**2*V) - Z*Rc**2/2)

! We prepare the values of the radial Fourier transform on a 3D grid
Vk(0) = 0
do i = 1, 3*(Ng/2+1)**2
    w = sqrt(real(i, dp))*dk
    Vk(i) = 4*pi/(L**3 * w**2) * (w*integrate(Rp, R*sin(w*R)*V) + Z*cos(w*Rc))
end do

! We fill out the 3D grid using the values from Vk
do i = 1, Ng
do j = 1, Ng
do k = 1, Ng
    idx = (i-1-Ng*nint((i-1.5_dp)/Ng))**2 + (j-1-Ng*nint((j-1.5_dp)/Ng))**2 &
            + (k-1-Ng*nint((k-1.5_dp)/Ng))**2
    VenG(i, j, k) = Vk(idx)
end do
end do
end do
end subroutine

subroutine radial_density_fourier(R, ne, L, neG)
! Takes radial density (given at origin) defined by:
!   ne(R) on a grid R
!   0 for r > maxval(R)
! and calculates a 3D Fourier transform of it into the neG grid. 'L' is the
! length of the box.
real(dp), intent(in) :: R(:), ne(:), L
real(dp), intent(out) :: neG(:, :, :)
integer :: Ng, i, j, k, idx
real(dp) :: Rp(size(R)), dk, Rc, w, nk(0:3*(size(neG, 1)/2+1)**2)
! Rp is the derivative of the mesh R'(t), which for uniform mesh is equal to
! the mesh step (rmax-rmin)/N:
Rp = (R(size(R)) - R(1)) / (size(R)-1)
Rc = R(size(R))
Ng = size(neG, 1)
call assert(size(neG, 2) == Ng)
call assert(size(neG, 3) == Ng)
dk = 2*pi/L

! We prepare the values of the radial Fourier transform on a 3D grid
nk(0) = 0
do i = 1, 3*(Ng/2+1)**2
    w = sqrt(real(i, dp))*dk
    nk(i) = 4*pi/(L**3 * w**2) * w*integrate(Rp, R*sin(w*R)*ne)
end do

! We fill out the 3D grid using the values from nk
do i = 1, Ng
do j = 1, Ng
do k = 1, Ng
    idx = (i-1-Ng*nint((i-1.5_dp)/Ng))**2 + (j-1-Ng*nint((j-1.5_dp)/Ng))**2 &
            + (k-1-Ng*nint((k-1.5_dp)/Ng))**2
    neG(i, j, k) = nk(idx)
end do
end do
end do
end subroutine

real(dp) pure function integrate(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)
s = integrate_trapz_1(Rp, f)
end function

subroutine free_energy_min(Nelec, Natom, L, G2, T_au, VenG, ne, energy_eps, &
        Eee, Een, Ts, Exc, Etot, cg_iter)
! Minimize the electronic free energy using the initial condition 'ne'. Returns
! the ground state in 'ne'. The free energy is returned in Etot, and it's
! components are returned in Eee, Een, Ts and Exc. The relation is:
!
!     Etot = Ts + Een + Eee + Exc
!
real(dp), intent(in) :: Nelec ! Number of electrons
integer, intent(in) :: Natom ! Number of atoms
real(dp), intent(in) :: L, G2(:, :, :), T_au, energy_eps
real(dp), intent(inout) :: ne(:, :, :)
complex(dp), intent(in) :: VenG(:, :, :)
real(dp), intent(out) :: Eee, Een, Ts, Exc, Etot
integer, intent(out) :: cg_iter ! # of CG iterations needed to converge

integer :: Ng
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

brent_eps = 1e-3_dp
max_iter = 2000
update_type = update_polak_ribiere
func_use_fft = .true.

last2 = 0
last3 = 0

Ng = size(ne, 1)

allocate(Hpsi(Ng, Ng, Ng))
allocate(psi(Ng, Ng, Ng))
allocate(psi_(Ng, Ng, Ng))
allocate(psi_prev(Ng, Ng, Ng))
allocate(phi(Ng, Ng, Ng))
allocate(phi_prime(Ng, Ng, Ng))
allocate(ksi(Ng, Ng, Ng))
allocate(ksi_prev(Ng, Ng, Ng))
allocate(eta(Ng, Ng, Ng))

psi = sqrt(ne)
psi_norm = integral(L, psi**2)
print *, "Initial norm of psi:", psi_norm
psi = sqrt(Nelec / psi_norm) * psi
psi_norm = integral(L, psi**2)
print *, "norm of psi:", psi_norm
! This returns H[n] = delta F / delta n, we save it to the Hpsi variable to
! save space:
call free_energy(L, G2, T_au, VenG, psi**2, Eee, Een, Ts, Exc, free_energy_, &
    Hpsi, calc_value=.false., calc_derivative=.true.)
! Hpsi = H[psi] = delta F / delta psi = 2*H[n]*psi, due to d/dpsi = 2 psi d/dn
Hpsi = Hpsi * 2*psi
mu = 1._dp / Nelec * integral(L, 0.5_dp * psi * Hpsi)
ksi = 2*mu*psi - Hpsi
phi = ksi
phi_prime = phi - 1._dp / Nelec *  integral(L, phi * psi) * psi
eta = sqrt(Nelec / integral(L, phi_prime**2)) * phi_prime
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
            call fourier2real(VenG, Ven)
            call precalc_C1C2C2(G2, psi, eta, C1, C2, C3)
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
    call free_energy(L, G2, T_au, VenG, psi**2, Eee, Een, Ts, Exc, &
        free_energy_, Hpsi, calc_value=.true., calc_derivative=.true.)
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
    print "('# ', i3, ' Etot/atom = ', f18.8, ' eV; last2 = ', es10.2, ' last3 = ',es10.2)", &
        iter, free_energy_ * Ha2eV / Natom, last2 * Ha2eV, last3 * Ha2eV
    if (iter > 3) then
        if (last3 < energy_eps) then
            ne = psi**2
            Etot = free_energy_
            cg_iter = iter
            return
        end if
    end if
    Hpsi = Hpsi * 2*psi ! d/dpsi = 2 psi d/dn
    mu = 1._dp / Nelec * integral(L, 0.5_dp * psi * Hpsi)
    ksi_prev = ksi
    ksi = 2*mu*psi - Hpsi
    select case(update_type)
        case(update_fletcher_reeves) ! Fletcher-Reeves
            gamma_n = integral(L, ksi**2)
        case(update_polak_ribiere)   ! Polak-Ribiere
            gamma_n = max(integral(L, ksi*(ksi-ksi_prev)), 0._dp)
        case default
            call stop_error("Unknown update type.")
    end select
    gamma_d = integral(L, ksi_prev**2)
    phi = ksi + gamma_n / gamma_d * phi
    phi_prime = phi - 1._dp / Nelec * integral(L, phi * psi) * psi
    eta = sqrt(Nelec / integral(L, phi_prime**2)) * phi_prime
    t2 = clock()
    print "('time: ', f6.3, ' fft: ', f6.3, ' (', f5.2, '%) # fft: ', i3, ' line search', f6.3)", &
        t2-t1, fft_time, fft_time / (t2-t1) * 100, fft_counter, t12-t11
end do
call stop_error("free_energy_minimization: The maximum number of iterations exceeded.")

contains

    real(dp) function func_nofft(theta) result(energy)
    real(dp), intent(in) :: theta
    call calc_F(L, T_au, Ven, C1, C2, C3, psi, eta, theta, energy)
    end function

    real(dp) function func_fft(theta) result(energy)
    real(dp), intent(in) :: theta
    psi_ = cos(theta) * psi + sin(theta) * eta
    call free_energy(L, G2, T_au, VenG, psi_**2, Eee, Een, Ts, Exc, &
        energy, Hpsi, calc_value=.true., calc_derivative=.false.)
    end function

end subroutine

subroutine forces(L, G, nG, Ven0G)
real(dp), intent(in) :: L, G(:, :, :, :)
complex(dp), intent(in) :: nG(:, :, :), Ven0G(:, :, :)
! Ven0G .... the potential returned from:
! radial_density_fourier(R, V, L, Z, Ven0G)
! x, y, z ... the ion position
real(dp) :: x, y, z, Fx, Fy, Fz
complex(dp) :: tmp(size(Ven0G, 1), size(Ven0G, 2), size(Ven0G, 3))
x = 1
y = 0
z = 0
tmp = i_ * nG * Ven0G * exp(-i_*(G(:,:,:,1)*x+G(:,:,:,2)*y+G(:,:,:,3)*z))
! Force: only take the real component of this:
Fx = real(sum(G(:,:,:,1) * L**3 * tmp), dp)
Fy = real(sum(G(:,:,:,2) * L**3 * tmp), dp)
Fz = real(sum(G(:,:,:,3) * L**3 * tmp), dp)
end subroutine

subroutine vtk_save(filename, Xn, f)
character(len=*), intent(in) :: filename
real(dp), intent(in) :: Xn(:, :, :, :), f(:, :, :)

integer :: Nx, Ny, Nz
integer :: u

Nx = size(Xn, 1)
Ny = size(Xn, 2)
Nz = size(Xn, 3)
call assert(size(Xn, 4) == 3)

call assert(size(f, 1) == Nx)
call assert(size(f, 2) == Ny)
call assert(size(f, 3) == Nz)

open(newunit=u, file=filename, status="replace")
write(u, "(a)") "# vtk DataFile Version 2.0"
write(u, "(a)") "Function at a real space uniform grid"
write(u, "(a)") "ASCII"
write(u, "(a)") "DATASET RECTILINEAR_GRID"
write(u, "(a, i8, i8, i8)") "DIMENSIONS", Nx, Ny, Nz

write(u, "(a, ' ', i8, ' ', a)") "X_COORDINATES", Nx, "float"
write(u, *) Xn(:, 1, 1, 1)
write(u, "(a, ' ', i8, ' ', a)") "Y_COORDINATES", Ny, "float"
write(u, *) Xn(1, :, 1, 2)
write(u, "(a, ' ', i8, ' ', a)") "Z_COORDINATES", Nz, "float"
write(u, *) Xn(1, 1, :, 3)

write(u, "(a, ' ', i10)") "POINT_DATA", size(f)
write(u, "(a)") "SCALARS MyFunction float"
write(u, "(a)") "LOOKUP_TABLE default"
write(u, *) f

close(u)
end subroutine

end module
