program test_free_energy_fft2

! Test for 4 Gaussian charges.

! This test uses FFT and produces the same result as test_free_energy4

use types, only: dp
use constants, only: i_
use ofdft, only: read_pseudo
use ofdft_fft, only: free_energy, radial_potential_fourier, &
    reciprocal_space_vectors, free_energy_min
use constants, only: Ha2eV, pi
use utils, only: loadtxt, stop_error, assert, linspace
use splines, only: spline3pars, iixmin, poly3, spline3ders
use interp3d, only: trilinear
use md, only: positions_fcc
implicit none
real(dp) :: Eee, Een, Ts, Exc, Etot
integer :: Ng
real(dp) :: Z
real(dp), allocatable :: R(:), G(:, :, :, :), G2(:, :, :)
real(dp), allocatable :: ne(:, :, :), dFdn(:, :, :)
real(dp), allocatable :: Ven0G(:, :, :)
complex(dp), allocatable :: VenG(:, :, :)
real(dp) :: L, T_eV, T_au
integer :: i, j, k
integer, parameter :: natom = 4
real(dp) :: X(3, natom), r_, alpha_ne, alpha_nen, x_, y_, z_

Ng = 80

L = 2
T_eV = 0.0862_dp
T_au = T_ev / Ha2eV

alpha_nen = 6
alpha_ne = 5

Z = 1

allocate(Ven0G(Ng, Ng, Ng), VenG(Ng, Ng, Ng), ne(Ng, Ng, Ng), dFdn(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng))
allocate(R(40000))
R = linspace(1._dp/40000, 0.9_dp, 40000)
call radial_potential_fourier(R, Z*erf(alpha_nen*R)/R, L, Z, Ven0G)

call reciprocal_space_vectors(L, G, G2)
ne = 1
do i = 1, Ng
do j = 1, Ng
do k = 1, Ng
    x_ = -L/2 + (i-1) * L / Ng
    y_ = -L/2 + (j-1) * L / Ng
    z_ = -L/2 + (k-1) * L / Ng
    r_ = sqrt(x_**2+y_**2+z_**2)
    ne(i, j, k) = 4*alpha_ne**3/pi**(3._dp/2) * exp(-alpha_ne**2*r_**2)
end do
end do
end do
call positions_fcc(X, L)
VenG = 0
do i = 1, natom
    VenG = VenG - Ven0G * exp(-i_ * &
        (G(:,:,:,1)*X(1,i) + G(:,:,:,2)*X(2,i) + G(:,:,:,3)*X(3,i)))
end do
call free_energy(L, G2, T_au, VenG, ne, Eee, Een, Ts, Exc, Etot, dFdn)
print *, "Summary of energies [a.u.]:"
print "('    Ts   = ', f14.8)", Ts
print "('    Een  = ', f14.8)", Een
print "('    Eee  = ', f14.8)", Eee
print "('    Exc  = ', f14.8)", Exc
print *, "   ---------------------"
print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV
print *, "Errors:"
print *, abs(Ts - (10.61904507_dp))
print *, abs(Een - (-2.92172113_dp))
print *, abs(Eee - (1.30109500_dp))
print *, abs(Exc - (-1.43805890_dp))
print *, abs(Etot - (7.56036004_dp))
call assert(abs(Ts - (10.61904507_dp)) < 1e-8_dp)
call assert(abs(Een - (-2.92172113_dp)) < 1e-8_dp)
call assert(abs(Eee - (1.30109500_dp)) < 1e-8_dp)
call assert(abs(Exc - (-1.43805890_dp)) < 1e-8_dp)
call assert(abs(Etot - (7.56036004_dp)) < 1e-8_dp)
end program
