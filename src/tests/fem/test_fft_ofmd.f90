program test_fft_ofmd
use types, only: dp
use constants, only: i_
use ofdft, only: read_pseudo
use ofdft_fft, only: free_energy, radial_density_fourier, &
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
real(dp), allocatable :: VenG0(:, :, :)
complex(dp), allocatable :: VenG(:, :, :)
real(dp) :: Rcut, L, T_eV, T_au
integer :: i, j, k
real(dp) :: x, y, z_, r_, alpha

Ng = 32

L = 2.997672536043746_dp
T_eV = 0.0862_dp
T_au = T_ev / Ha2eV

call read_pseudo("H.pseudo.gaussian2", R, V, Z, Ediff)
Rcut = R(size(R))

allocate(VenG0(Ng, Ng, Ng), VenG(Ng, Ng, Ng), ne(Ng, Ng, Ng), dFdn(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng))
call radial_density_fourier(R, V, L, Z, VenG0)

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
VenG = -VenG0 * exp(i_*(G(:,:,:,1)*x+G(:,:,:,2)*y+G(:,:,:,3)*z_))
call free_energy(L, G2, T_au, VenG, ne, Eee, Een, Ts, Exc, Etot, dFdn)
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
call free_energy_min(1, L, G2, T_au, VenG, ne, Eee, Een, Ts, Exc, Etot)
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
