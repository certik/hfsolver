program test_free_energy_fft

! nuclear charge: 1 Gaussian
! electronic charge: 1 Gaussian
! calculation: single free energy evaluation

! This test uses FFT and produces the same result as test_free_energy3

use types, only: dp
use constants, only: i_
use ofdft, only: read_pseudo
use ofdft_fft, only: free_energy, radial_potential_fourier, &
    reciprocal_space_vectors, free_energy_min
use constants, only: Ha2eV, pi
use utils, only: loadtxt, stop_error, assert, linspace
use splines, only: spline3pars, iixmin, poly3, spline3ders
use interp3d, only: trilinear
use converged_energies, only: one_gaussian
implicit none
real(dp) :: Eee, Een, Ts, Exc, Etot, Etot_conv
integer :: Ng
real(dp) :: Z
real(dp), allocatable :: R(:), G(:, :, :, :), G2(:, :, :)
real(dp), allocatable :: ne(:, :, :), dFdn(:, :, :)
real(dp), allocatable :: VenG0(:, :, :)
real(dp) :: V0
complex(dp), allocatable :: VenG(:, :, :)
real(dp) :: L, T_eV, T_au
integer :: i, j, k
real(dp) :: x, y, z_, r_, alpha_ne, alpha_nen

Ng = 80

L = 2
T_eV = 0.0862_dp
T_au = T_ev / Ha2eV

alpha_nen = 6
alpha_ne = 5

Z = 1

allocate(VenG0(Ng, Ng, Ng), VenG(Ng, Ng, Ng), ne(Ng, Ng, Ng), dFdn(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng))
allocate(R(10000))
R = linspace(1._dp/10000, 0.9_dp, 10000)
call radial_potential_fourier(R, Z*erf(alpha_nen*R)/R, L, Z, VenG0, V0)

call reciprocal_space_vectors(L, G, G2)
ne = 1
do i = 1, Ng
do j = 1, Ng
do k = 1, Ng
    x = -L/2 + (i-1) * L / Ng
    y = -L/2 + (j-1) * L / Ng
    z_ = -L/2 + (k-1) * L / Ng
    r_ = sqrt(x**2+y**2+z_**2)
    ne(i, j, k) = alpha_ne**3/pi**(3._dp/2) * exp(-alpha_ne**2*r_**2)
end do
end do
end do
! Ion position:
x = L/2 + L/64
y = L/2
z_ = L/2
VenG = -VenG0 * exp(-i_*(G(:,:,:,1)*x+G(:,:,:,2)*y+G(:,:,:,3)*z_))
call free_energy(L, G2, T_au, VenG, ne, Eee, Een, Ts, Exc, Etot, dFdn)
Etot_conv = sum(one_gaussian)
print *, "Summary of energies [a.u.]:"
print "('    Ts   = ', f14.8)", Ts
print "('    Een  = ', f14.8)", Een
print "('    Eee  = ', f14.8)", Eee
print "('    Exc  = ', f14.8)", Exc
print *, "   ---------------------"
print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV
print *, "Errors:"
print *, abs(Ts - one_gaussian(1))
print *, abs(Een - one_gaussian(2))
print *, abs(Eee - one_gaussian(3))
print *, abs(Exc - one_gaussian(4))
print *, abs(Etot - Etot_conv)
call assert(abs(Ts - one_gaussian(1)) < 1e-8_dp)
call assert(abs(Een - one_gaussian(2)) < 1e-8_dp)
call assert(abs(Eee - one_gaussian(3)) < 1e-8_dp)
call assert(abs(Exc - one_gaussian(4)) < 1e-8_dp)
call assert(abs(Etot - Etot_conv) < 1e-8_dp)
end program
