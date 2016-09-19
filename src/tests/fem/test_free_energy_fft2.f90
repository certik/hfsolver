program test_free_energy_fft2

! nuclear charge: 1 Gaussian
! electronic charge: 4 Gaussians
! calculation: single free energy evaluation

! This test uses FFT and produces the same result as test_free_energy4

use types, only: dp
use constants, only: i_
use ofdft, only: read_pseudo
use ofdft_fft, only: free_energy, radial_potential_fourier, &
    reciprocal_space_vectors, free_energy_min, real2fourier, real_space_vectors
use constants, only: Ha2eV, pi
use utils, only: loadtxt, stop_error, assert, linspace
use splines, only: spline3pars, iixmin, poly3, spline3ders
use interp3d, only: trilinear
use md, only: positions_fcc
use converged_energies, only: four_gaussians, four_gaussians_min
implicit none
real(dp) :: Eee, Een, Ts, Exc, Etot, Etot_conv
integer :: Ng
real(dp) :: Z
real(dp), allocatable :: R(:), G(:, :, :, :), G2(:, :, :), Xn(:, :, :, :)
real(dp), allocatable :: ne(:, :, :), dFdn(:, :, :)
real(dp), allocatable :: Ven0G(:, :, :), fac(:, :, :)
real(dp) :: V0
complex(dp), allocatable :: VenG(:, :, :), neG(:, :, :)
real(dp) :: L, T_eV, T_au
integer :: i, j, k
integer, parameter :: natom = 4
real(dp) :: X(3, natom), fen(3, natom), r_, alpha_ne, alpha_nen, x_, y_, z_
integer :: cg_iter

Ng = 80

L = 2
T_eV = 0.0862_dp
T_au = T_ev / Ha2eV

alpha_nen = 6
alpha_ne = 5

Z = 1

allocate(Ven0G(Ng, Ng, Ng), VenG(Ng, Ng, Ng), ne(Ng, Ng, Ng), dFdn(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng), fac(Ng, Ng, Ng), neG(Ng, Ng, Ng))
allocate(R(40000), Xn(Ng, Ng, Ng, 3))
R = linspace(1._dp/40000, 0.9_dp, 40000)
call radial_potential_fourier(R, Z*erf(alpha_nen*R)/R, L, Z, Ven0G, V0)

call real_space_vectors([L, L, L], Xn)
call reciprocal_space_vectors([L, L, L], G, G2)
ne = 1
do i = 1, Ng
do j = 1, Ng
do k = 1, Ng
    x_ = Xn(i, j, k, 1) - L/2
    y_ = Xn(i, j, k, 2) - L/2
    z_ = Xn(i, j, k, 3) - L/2
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
call free_energy(L, G2, T_au, VenG, ne, Eee, Een, Ts, Exc, Etot, dFdn, &
    calc_value=.true., calc_derivative=.false.)
Etot_conv = sum(four_gaussians)
print *, "Summary of energies [a.u.]:"
print "('    Ts   = ', f14.8)", Ts
print "('    Een  = ', f14.8)", Een
print "('    Eee  = ', f14.8)", Eee
print "('    Exc  = ', f14.8)", Exc
print *, "   ---------------------"
print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV
print *, "Errors:"
print *, abs(Ts - four_gaussians(1))
print *, abs(Een - four_gaussians(2))
print *, abs(Eee - four_gaussians(3))
print *, abs(Exc - four_gaussians(4))
print *, abs(Etot - Etot_conv)
call assert(abs(Ts - four_gaussians(1)) < 1e-8_dp)
call assert(abs(Een - four_gaussians(2)) < 1e-8_dp)
call assert(abs(Eee - four_gaussians(3)) < 1e-8_dp)
call assert(abs(Exc - four_gaussians(4)) < 1e-8_dp)
call assert(abs(Etot - Etot_conv) < 1e-8_dp)

! Minimization

ne = 4 / L**3
call free_energy_min(4._dp, 4, L, G2, T_au, VenG, ne, 1e-12_dp, Eee, Een, Ts, &
    Exc, Etot, cg_iter)

Etot_conv = sum(four_gaussians_min)
print *, "Summary of energies [a.u.]:"
print "('    Ts   = ', f14.8)", Ts
print "('    Een  = ', f14.8)", Een
print "('    Eee  = ', f14.8)", Eee
print "('    Exc  = ', f14.8)", Exc
print *, "   ---------------------"
print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV
print *, "Errors:"
print *, abs(Ts - four_gaussians_min(1))
print *, abs(Een - four_gaussians_min(2))
print *, abs(Eee - four_gaussians_min(3))
print *, abs(Exc - four_gaussians_min(4))
print *, abs(Etot - Etot_conv)
call assert(abs(Ts - four_gaussians_min(1)) < 5e-7_dp)
call assert(abs(Een - four_gaussians_min(2)) < 5e-7_dp)
call assert(abs(Eee - four_gaussians_min(3)) < 5e-8_dp)
call assert(abs(Exc - four_gaussians_min(4)) < 5e-8_dp)
call assert(abs(Etot - Etot_conv) < 1e-10_dp)

! ----------------------------------------------------------------------
! Forces calculation:

call real2fourier(ne, neG)
fen = 0
do i = 1, natom
    fac = L**3*Ven0G*aimag(neG*exp(-i_ * &
        (G(:,:,:,1)*X(1,i) + G(:,:,:,2)*X(2,i) + G(:,:,:,3)*X(3,i))))
    fen(1, i) = sum(G(:,:,:,1)*fac)
    fen(2, i) = sum(G(:,:,:,2)*fac)
    fen(3, i) = sum(G(:,:,:,3)*fac)
end do

print *, "forces FFT:"
print *, fen(:, 1)
print *, fen(:, 2)
print *, fen(:, 3)
print *, fen(:, 4)
end program
