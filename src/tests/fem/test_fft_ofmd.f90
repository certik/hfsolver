module test_fft_ofmd_utils

! OFMD using FFT. We use 3D reciprocal grid (:, :, :). Variables in reciprocal
! (Fourier) space have F appended to them. E.g. Ven is the electron-nucleus
! potential in real space, VenF is the same potential in reciprocal space.

use types, only: dp
use utils, only: assert, zeros, stop_error
use constants, only: pi
use xc, only: xc_pz
use optimize, only: bracket, brent
use fourier, only: fft3_inplace, ifft3_inplace
implicit none
private
public read_pseudo, radial_density_fourier, free_energy, &
    reciprocal_space_vectors, free_energy_min

! Update types for nonlinear conjugate gradient method:
integer, parameter :: update_fletcher_reeves = 1
integer, parameter :: update_polak_ribiere   = 2

contains

subroutine reciprocal_space_vectors(L, G, G2)
real(dp), intent(in) :: L
! G(:, :, :, i) where i=1, 2, 3 are the x, y, z components
! G2(:, :, :) are the squares of G
real(dp), intent(out) :: G(:, :, :, :), G2(:, :, :)
real(dp) :: G1D(size(G, 1))
integer :: Ng, i, j, k
Ng = size(G, 1)
forall(i=1:Ng) G1D(i) = 2*pi/L * (i-1-Ng*nint((i-1.5_dp)/Ng))
forall(i=1:Ng, j=1:Ng, k=1:Ng)
    G(i, j, k, 1) = G1D(i)
    G(i, j, k, 2) = G1D(j)
    G(i, j, k, 3) = G1D(k)
end forall
G(1, 1, 1, :) = 1 ! To avoid division by 0
G2 = G(:,:,:,1)**2 + G(:,:,:,2)**2 + G(:,:,:,3)**2
end subroutine

subroutine real2fourier(x, xF)
! Calculates Discrete Fourier Transform, with the same normalization as the
! Fourier Transform for periodic systems (which is Fourier Series).
real(dp), intent(in) :: x(:, :, :)
complex(dp), intent(out) :: xF(:, :, :)
integer :: Ng
Ng = size(x, 1)
xF = x
call fft3_inplace(xF) ! Calculates sum_{n=0}^{N-1} e^{-2*pi*i*k*n/N}*x(n)
xF = xF / size(x)     ! The proper normalization is to divide by N
end subroutine

subroutine fourier2real(xF, x)
! Calculates Inverse Discrete Fourier Transform. xF must follow the same
! normalization as defined by real2fourier(), i.e. in the following calls, 'x2'
! will be equal to 'x':
! call real2fourier(x, xF)
! call fourier2real(xF, x2)
complex(dp), intent(in) :: xF(:, :, :)
real(dp), intent(out) :: x(:, :, :)
complex(dp) :: tmp(size(xF,1), size(xF,2), size(xF,3))
integer :: Ng
Ng = size(x, 1)
tmp = xF
call ifft3_inplace(tmp) ! Calculates sum_{k=0}^{N-1} e^{2*pi*i*k*n/N}*X(k)
x = real(tmp, dp)       ! The result is already normalized
! TODO: make this strict:
call assert(maxval(abs(aimag(tmp))) < 1e-1_dp)
end subroutine

real(dp) function integral(L, f) result(r)
real(dp), intent(in) :: L, f(:, :, :)
r = sum(f) * L**3 / size(f)
end function

real(dp) function integralF(L, fF) result(r)
real(dp), intent(in) :: L
complex(dp), intent(in) :: fF(:, :, :)
complex(dp) :: s
s = sum(fF) * L**3
r = real(s, dp)
call assert(abs(aimag(s)) < 1e-12_dp)
end function

real(dp) function integralF2(fF, L) result(r)
real(dp), intent(in) :: L, fF(:, :, :)
r = sum(fF) * L**3
end function

subroutine parabola_vertex(x1, y1, x2, y2, x3, y3, xv, yv)
real(dp), intent(in) :: x1, y1, x2, y2, x3, y3
real(dp), intent(out) :: xv, yv
real(dp) :: denom, A, B, C
denom = (x1 - x2) * (x1 - x3) * (x2 - x3)
A     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
B     = (x3**2 * (y1 - y2) + x2**2 * (y3 - y1) + x1**2 * (y2 - y3)) / denom
C     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + &
            x1 * x2 * (x1 - x2) * y3) / denom
xv = -B / (2*A)
yv = C - B**2 / (4*A)
end subroutine

subroutine free_energy(L, G2, T_au, VenF, ne, Eee, Een, Ts, Exc, Etot, dFdn)
real(dp), intent(in) :: L, G2(:, :, :), T_au, ne(:, :, :)
complex(dp), intent(in) :: VenF(:, :, :)
real(dp), intent(out) :: Eee, Een, Ts, Exc, Etot

! dFdn returns "delta F / delta n", the functional derivative with respect to
! the density "n". Use the relation
!     d/dpsi = 2 psi d/dn
! to obtain the derivative with respect to psi (i.e. multiply Hpsi by 2*psi).
real(dp), intent(out) :: dFdn(:, :, :)

real(dp), dimension(size(VenF,1), size(VenF,2), size(VenF,3)) :: y, F0, &
    exc_density, Vee, Ven, Vxc, dF0dn
complex(dp), dimension(size(VenF,1), size(VenF,2), size(VenF,3)) :: &
    neF, VeeF
real(dp) :: beta, dydn
integer :: i, j, k, Ng
Ng = size(VenF, 1)

call real2fourier(ne, neF)

VeeF = 4*pi*neF / G2
VeeF(1, 1, 1) = 0

! Hartree energy
!Eee = integralF2(L, real(VeeF)*real(neF)+aimag(VeeF)*aimag(neF)) / 2
Eee = integralF(L, VeeF*conjg(neF)) / 2
! Electron-nucleus energy
!Een = integralF2(L, real(VenF)*real(neF)+aimag(VenF)*aimag(neF))
Een = integralF(L, VenF*conjg(neF))

! Kinetic energy using Perrot parametrization
beta = 1/T_au
! The density must be positive, the f(y) fails for negative "y". Thus we use
! ne.
y = pi**2 / sqrt(2._dp) * beta**(3._dp/2) * ne
if (any(y < 0)) call stop_error("Density must be positive")
F0 = ne / beta * f(y)
Ts = integral(L, F0)
! Exchange and correlation potential
do k = 1, Ng
do j = 1, Ng
do i = 1, Ng
    call xc_pz(ne(i, j, k), exc_density(i, j, k), Vxc(i, j, k))
end do
end do
end do
Exc = integral(L, exc_density * ne)
Etot = Ts + Een + Eee + Exc

! Calculate the derivative
dydn = pi**2 / sqrt(2._dp) * beta**(3._dp/2)
! F0 = ne / beta * f(y)
! d F0 / d n =
dF0dn = 1 / beta * f(y) + ne / beta * f(y, deriv=.true.) * dydn

call fourier2real(VenF, Ven)
call fourier2real(VeeF, Vee)

!print *, dF0dn
!print *, Vee
!print *, Ven
!print *, Vxc
!print *, "---------"
dFdn = dF0dn + Vee + Ven + Vxc
dFdn = dFdn * L**3
end subroutine

real(dp) elemental function f(y, deriv)
! Function f(y) from Appendix A in [1].
!
! [1] Perrot, F. (1979). Gradient correction to the statistical electronic free
! energy at nonzero temperatures: Application to equation-of-state
! calculations. Physical Review A, 20(2), 586–594.
real(dp), intent(in) :: y ! must be positive
! if deriv == .true. compute df/dy instead. Default .false.
logical, intent(in), optional :: deriv
real(dp), parameter :: y0 = 3*pi/(4*sqrt(2._dp))
real(dp), parameter :: c(*) = [-0.8791880215_dp, 0.1989718742_dp, &
    0.1068697043e-2_dp, -0.8812685726e-2_dp, 0.1272183027e-1_dp, &
    -0.9772758583e-2_dp, 0.3820630477e-2_dp, -0.5971217041e-3_dp]
real(dp), parameter :: d(*) = [0.7862224183_dp, -0.1882979454e1_dp, &
    0.5321952681_dp, 0.2304457955e1_dp, -0.1614280772e2_dp, &
    0.5228431386e2_dp, -0.9592645619e2_dp, 0.9462230172e2_dp, &
    -0.3893753937e2_dp]
real(dp) :: u
integer :: i
logical :: deriv_
deriv_ = .false.
if (present(deriv)) deriv_ = deriv

if (.not. deriv_) then
    if (y <= y0) then
        f = log(y)
        do i = 0, 7
            f = f + c(i+1) * y**i
        end do
    else
        u = y**(2._dp / 3)
        f = 0
        do i = 0, 8
            f = f + d(i+1) / u**(2*i-1)
        end do
        ! Note: Few terms in [1] have "y" instead of "u" in them for y > y0, but
        ! that is obviously a typo.
    end if
else
    if (y <= y0) then
        f = 1 / y
        do i = 0, 6
            f = f + (i+1) * c(i+2) * y**i
        end do
    else
        u = y**(2._dp / 3)
        f = 0
        do i = 0, 8
            f = f - (2*i-1) * d(i+1) / u**(2*i)
        end do
        f = f * 2._dp/3 / y**(1._dp/3)
    end if
end if
end function


subroutine read_pseudo(filename, R, V, Z, Ediff)
! Reads the pseudopotential from the file 'filename'.
character(len=*), intent(in) :: filename   ! File to read from, e.g. "H.pseudo"
real(dp), allocatable, intent(out) :: R(:) ! radial grid [0, Rcut]
! potential on the radial grid. The potential smoothly changes into -1/R for
! r > Rcut, where Rcut = R(size(R)) is the cut-off radius
real(dp), allocatable, intent(out) :: V(:)
real(dp), intent(out) :: Z     ! Nuclear charge
real(dp), intent(out) :: Ediff ! The energy correction
real(dp) :: Rcut
integer :: N, i, u
open(newunit=u, file=filename, status="old")
read(u, *) Z, N, Rcut, Ediff
allocate(R(N), V(N))
do i = 1, N
    read(u, *) R(i), V(i)
end do
close(u)
R = R * Rcut
end subroutine

subroutine radial_density_fourier(R, V, L, Z, VenF)
real(dp), intent(in) :: R(:), V(:), L, Z
complex(dp), intent(out) :: VenF(:, :, :)
integer :: Ng, i, j, k, idx
real(dp) :: Rp(size(R)), dk, Rc, w, Vk(0:3*(size(VenF, 1)/2+1)**2)
! Rp is the derivative of the mesh R'(t), which for uniform mesh is equal to
! the mesh step (rmax-rmin)/N:
Rp = (R(size(R)) - R(1)) / (size(R)-1)
Rc = R(size(R))
Ng = size(VenF, 1)
call assert(size(VenF, 2) == Ng)
call assert(size(VenF, 3) == Ng)
dk = 2*pi/L
! We prepare the values of the radial Fourier transform on a 3D grid
Vk(0) = 0
do i = 1, 3*(Ng/2+1)**2
    w = sqrt(real(i, dp))*dk
    Vk(i) = 4*pi/(L**3 * w**2) * (w*integrate(Rp, R*sin(w*R)*V) + cos(w*Rc))
end do

! We fill out the 3D grid using the values from Vk
do i = 1, Ng
do j = 1, Ng
do k = 1, Ng
    idx = (i-1-Ng*nint((i-1.5_dp)/Ng))**2 + (j-1-Ng*nint((j-1.5_dp)/Ng))**2 &
            + (k-1-Ng*nint((k-1.5_dp)/Ng))**2
    VenF(i, j, k) = Vk(idx) * Z
end do
end do
end do
end subroutine

real(dp) pure function integrate(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)
s = integrate_trapz_1(Rp, f)
end function

real(dp) pure function integrate_trapz_1(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)
real(dp) :: g(size(Rp))
integer :: N
N = size(Rp)
g = f * Rp
s = (g(1) + g(N)) / 2
s = s + sum(g(2:N-1))
end function

subroutine free_energy_min(L, G2, T_au, VenF, ne, Eee, Een, Ts, Exc, Etot)
real(dp), intent(in) :: L, G2(:, :, :), T_au
real(dp), intent(inout) :: ne(:, :, :)
complex(dp), intent(in) :: VenF(:, :, :)
real(dp), intent(out) :: Eee, Een, Ts, Exc, Etot

integer :: Ng
real(dp), allocatable :: free_energies(:)
real(dp), allocatable, dimension(:, :, :) :: Hpsi, &
    psi, psi_, psi_prev, ksi, ksi_prev, phi, phi_prime, eta
integer :: iter, max_iter
real(dp) :: mu, energy_eps, last3, brent_eps, free_energy_, &
    gamma_d, gamma_n, theta, theta_a, theta_b, theta_c, fa, fb, fc
real(dp) :: f2
real(dp) :: Nelec
real(dp) :: psi_norm
integer :: update_type
!real(dp) :: A, B

!energy_eps = 3.6749308286427368e-5_dp
energy_eps = 1e-9_dp
brent_eps = 1e-3_dp
max_iter = 200
update_type = update_polak_ribiere

Nelec = 1 ! One electron
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
call free_energy(L, G2, T_au, VenF, psi**2, Eee, Een, Ts, Exc, free_energy_, Hpsi)
! Hpsi = H[psi] = delta F / delta psi = 2*H[n]*psi, due to d/dpsi = 2 psi d/dn
Hpsi = Hpsi * 2*psi
mu = 1._dp / Nelec * integral(L, 0.5_dp * psi * Hpsi)
ksi = 2*mu*psi - Hpsi
phi = ksi
phi_prime = phi - 1._dp / Nelec *  integral(L, phi * psi) * psi
eta = sqrt(Nelec / integral(L, phi_prime**2)) * phi_prime
theta = pi/2
print *, "Summary of energies [a.u.]:"
print "('    Ts   = ', f14.8)", Ts
print "('    Een  = ', f14.8)", Een
print "('    Eee  = ', f14.8)", Eee
print "('    Exc  = ', f14.8)", Exc
print *, "   ---------------------"
print "('    Etot = ', f14.8, ' a.u.')", free_energy_
allocate(free_energies(max_iter))
gamma_n = 0
do iter = 1, max_iter
    !Hpsi = Hpsi/(2*psi)
    ! Formula (36) in Jiang & Yang
    !A = integral(L, psi*Hpsi*psi) - integral(L, eta*Hpsi*eta)
    !B = 2*integral(L, eta*Hpsi*psi)
    !print *, "theta? =", 0.5_dp * atan(B/A)
    theta_a = 0
    theta_b = mod(theta, 2*pi)
    call bracket(func, theta_a, theta_b, theta_c, fa, fb, fc, 100._dp, 20, verbose=.false.)
    if (iter < 2) then
        call brent(func, theta_a, theta_b, theta_c, brent_eps, 50, theta, &
            free_energy_, verbose=.true.)
    else
        call parabola_vertex(theta_a, fa, theta_b, fb, theta_c, fc, theta, f2)
    end if
    ! TODO: We probably don't need to recalculate free_energy_ here:
    psi_prev = psi
    psi = cos(theta) * psi + sin(theta) * eta
    call free_energy(L, G2, T_au, VenF, psi**2, Eee, Een, Ts, Exc, &
        free_energy_, Hpsi)
    print *, "Iteration:", iter
    psi_norm = integral(L, psi**2)
    print *, "Norm of psi:", psi_norm
    print *, "mu =", mu
    print *, "|ksi| =", sqrt(gamma_n)
    print *, "theta =", theta
    print *, "Summary of energies [a.u.]:"
    print "('    Ts   = ', f14.8)", Ts
    print "('    Een  = ', f14.8)", Een
    print "('    Eee  = ', f14.8)", Eee
    print "('    Exc  = ', f14.8)", Exc
    print *, "   ---------------------"
    print "('    Etot = ', f14.8, ' a.u.')", free_energy_
    free_energies(iter) = free_energy_
    if (iter > 3) then
        last3 = maxval(free_energies(iter-3:iter)) - &
            minval(free_energies(iter-3:iter))
        if (last3 < energy_eps) then
            ne = psi**2
            Etot = free_energy_
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
end do
call stop_error("free_energy_minimization: The maximum number of iterations exceeded.")

contains

    real(dp) function func(theta) result(energy)
    real(dp), intent(in) :: theta
    energy = theta
    psi_ = cos(theta) * psi + sin(theta) * eta
    call free_energy(L, G2, T_au, VenF, psi_**2, Eee, Een, Ts, Exc, &
        energy, Hpsi)
    end function

end subroutine


end module


! ------------------------------------------------------------------------


program test_fft_ofmd
use types, only: dp
use constants, only: i_
use test_fft_ofmd_utils, only: free_energy, read_pseudo, &
    radial_density_fourier, &
    reciprocal_space_vectors, free_energy_min
use constants, only: Ha2eV, pi
use utils, only: loadtxt, stop_error, assert
use splines, only: spline3pars, iixmin, poly3, spline3ders
use interp3d, only: trilinear
implicit none
real(dp) :: Eee, Een, Ts, Exc, Etot
integer :: Ng
real(dp) :: Z, Ediff
real(dp), allocatable :: R(:), V(:), G(:, :, :, :), G2(:, :, :)
real(dp), allocatable :: ne(:, :, :), dFdn(:, :, :)
complex(dp), allocatable :: VenF(:, :, :)
real(dp) :: Rcut, L, T_eV, T_au
integer :: i, j, k
real(dp) :: x, y, z_, r_, alpha

Ng = 32

L = 2.997672536043746_dp
T_eV = 0.0862_dp
T_au = T_ev / Ha2eV

call read_pseudo("H.pseudo.gaussian", R, V, Z, Ediff)
Rcut = R(size(R))

allocate(VenF(Ng, Ng, Ng), ne(Ng, Ng, Ng), dFdn(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng))
call radial_density_fourier(R, V, L, Z, VenF)

call reciprocal_space_vectors(L, G, G2)
ne = 1
do i = 1, Ng
do j = 1, Ng
do k = 1, Ng
    x = -L/2 + (i-1) * L / Ng
    y = -L/2 + (j-1) * L / Ng
    z_ = -L/2 + (k-1) * L / Ng
    r_ = sqrt(x**2+y**2+z_**2)
    alpha = 3
    ne(i, j, k) = alpha**3/pi**(3._dp/2) * exp(-alpha**2*r_**2)
end do
end do
end do
! Ion position:
x = L/2 + L/64
y = L/2
z_ = L/2
VenF = -VenF * exp(i_*(G(:,:,:,1)*x+G(:,:,:,2)*y+G(:,:,:,3)*z_))
call free_energy(L, G2, T_au, VenF, ne, Eee, Een, Ts, Exc, Etot, dFdn)
print *, "Summary of energies [a.u.]:"
print "('    Ts   = ', f14.8)", Ts
print "('    Een  = ', f14.8)", Een
print "('    Eee  = ', f14.8)", Eee
print "('    Exc  = ', f14.8)", Exc
print *, "   ---------------------"
print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV
call assert(abs(Ts - (3.82284833)) < 5e-8)
call assert(abs(Een - (-2.33127597)) < 5e-8)
call assert(abs(Eee - (0.73653527)) < 5e-8)
call assert(abs(Exc - (-0.88363737)) < 5e-8)
call assert(abs(Etot - (1.34447026)) < 1e-8)
ne=1._dp / L**3
call free_energy_min(L, G2, T_au, VenF, ne, Eee, Een, Ts, Exc, Etot)
print *, "Ng =", Ng
print *, "Rcut =", Rcut
print *, "T_au =", T_au
print *, "Summary of energies [a.u.]:"
print "('    Ts   = ', f14.8)", Ts
print "('    Een  = ', f14.8)", Een
print "('    Eee  = ', f14.8)", Eee
print "('    Exc  = ', f14.8)", Exc
print *, "   ---------------------"
print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV
call assert(abs(Ts - (0.75568820)) < 1e-4)
call assert(abs(Een - (-0.78524779)) < 1e-4)
call assert(abs(Eee - (0.07603492)) < 1e-4)
call assert(abs(Exc - (-0.38053818)) < 1e-4)
call assert(abs(Etot - (-0.33406285)) < 1e-8)
end program
