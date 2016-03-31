program test_functional_derivative_min_fft

! Tests the functional derivative by direct minimization (imaginary time)
! versus CG.

use types, only: dp
use constants, only: i_
use ofdft, only: read_pseudo
use ofdft_fft, only: free_energy, radial_potential_fourier, &
    reciprocal_space_vectors, free_energy_min, real2fourier, integral, &
    fourier2real, real_space_vectors
use constants, only: Ha2eV
use utils, only: loadtxt, stop_error, assert, linspace
use splines, only: spline3pars, iixmin, poly3, spline3ders
use interp3d, only: trilinear
use md, only: positions_fcc, positions_bcc
implicit none
real(dp) :: Eee, Een, Ts, Exc, Etot, Etot_conv, mu_conv
integer :: Ng
real(dp) :: Z
real(dp), allocatable :: R(:), G(:, :, :, :), G2(:, :, :), psi(:,:,:)
real(dp), allocatable :: ne(:, :, :), Hn(:, :, :), Hn2(:,:,:), Ven_rad(:)
real(dp), allocatable :: Ven0G(:, :, :), current(:,:,:,:), Xn(:,:,:,:), X(:,:)
real(dp) :: V0
complex(dp), allocatable, dimension(:,:,:) :: VenG, tmp
complex(dp), allocatable, dimension(:,:,:,:) :: dpsi
real(dp) :: L, T_eV, T_au
integer :: i, u
integer :: natom = 128
real(dp) :: alpha_nen, mu, mu_Hn, dt, psi_norm, Ediff
integer :: cg_iter
real(dp) :: E0, t, omega

Ng = 32

L = 8.1049178668765851_dp
T_eV = 34.5_dp
T_au = T_ev / Ha2eV

alpha_nen = 6

allocate(Ven0G(Ng, Ng, Ng), VenG(Ng, Ng, Ng), ne(Ng, Ng, Ng), Hn(Ng, Ng, Ng))
allocate(Hn2(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng), psi(Ng, Ng, Ng))
allocate(dpsi(Ng, Ng, Ng, 3), current(Ng, Ng, Ng, 3), tmp(Ng, Ng, Ng))
allocate(Xn(Ng, Ng, Ng, 3))
print *, "Load initial position"
allocate(X(3, natom))
call positions_bcc(X, L)
print *, "Radial nuclear potential FFT"
call read_pseudo("D.pseudo", R, Ven_rad, Z, Ediff)
call radial_potential_fourier(R, Ven_rad, L, Z, Ven0G, V0)
print *, "    Done."

call real_space_vectors(L, Xn)
call reciprocal_space_vectors(L, G, G2)
VenG = 0
do i = 1, natom
    VenG = VenG - Ven0G * exp(-i_ * &
        (G(:,:,:,1)*X(1,i) + G(:,:,:,2)*X(2,i) + G(:,:,:,3)*X(3,i)))
end do

! Minimization

ne = natom / L**3

call free_energy(L, G2, T_au, VenG, ne, Eee, Een, Ts, Exc, Etot, Hn, &
    calc_value=.true., calc_derivative=.true.)

mu = sum(Hn)/size(Hn)
print *, "mu = ", mu
print *, "max(abs(H-mu)) = ", maxval(abs(Hn - mu))
call assert(all(ne > 0))

print *
print *, "------------------------------------------------------------------"
print *, "Propagation"

! Propagate

print *, "E_max =", maxval(abs(Hn)), "; dt <", 1/maxval(abs(Hn))
dt = 1/maxval(abs(Hn)) / 10 ! set dt 10x smaller than the limit
dt = 1e-4_dp
print *, "dt =", dt

! Do first step by hand:
print *, "First step"
psi = sqrt(ne)

t = 0

psi = psi - dt*Hn*psi

ne = psi**2
psi_norm = integral(L, ne)
print *, "Initial norm of psi:", psi_norm
psi = sqrt(natom / psi_norm) * psi
ne = psi**2
psi_norm = integral(L, ne)
print *, "norm of psi:", psi_norm

E0 = 1e-3_dp
omega = 0.05

open(newunit=u, file="log.txt", status="replace")
close(u)

do i = 1, 200
    t = t + dt
    print *, "iter =", i, "time =", t
    psi = psi - dt*Hn*psi
    ne = psi**2
    psi_norm = integral(L, ne)
    print *, "Initial norm of psi:", psi_norm
    psi = sqrt(natom / psi_norm) * psi
    ne = psi**2
    psi_norm = integral(L, ne)
    print *, "norm of psi:", psi_norm

    call free_energy(L, G2, T_au, VenG, ne, Eee, Een, Ts, Exc, Etot, Hn, &
        calc_value=.true., calc_derivative=.true.)
    Etot = Ts + Een + Eee + Exc

    mu = 1._dp / natom * integral(L, ne * Hn)
    mu_Hn = sum(Hn)/size(Hn)
    print *, mu, mu_Hn
    print *, "Summary of energies [a.u.]:"
    print "('    Ts   = ', f14.8)", Ts
    print "('    Een  = ', f14.8)", Een
    print "('    Eee  = ', f14.8)", Eee
    print "('    Exc  = ', f14.8)", Exc
    print *, "   ---------------------"
    print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV


    open(newunit=u, file="log.txt", position="append", status="old")
    write(u, *) i, Etot, mu, mu_Hn
    close(u)

end do
print *, "Done"

Etot_conv = -172.12475770606159_dp
mu_conv = 96.415580964855209_dp

print *, "Etot:", Etot, abs(Etot - Etot_conv)
print *, "mu:", mu, abs(mu - mu_conv)
print *, "mu_Hn:", mu_Hn, abs(mu_Hn - mu_conv)
print *, "abs(mu-mu_Hn):", abs(mu - mu_Hn)
call assert(abs(Etot - Etot_conv) < 1e-12_dp)
call assert(abs(mu - mu_conv) < 1e-11_dp)
call assert(abs(mu - mu_Hn) < 1e-10_dp)

! Now compare against CG minimization
print *, "CG minimization:"
ne = natom / L**3
call free_energy_min(real(natom, dp), natom, L, G2, T_au, VenG, ne, &
    1e-12_dp, Eee, Een, Ts, Exc, Etot, cg_iter)
call free_energy(L, G2, T_au, VenG, ne, Eee, Een, Ts, Exc, Etot, Hn, &
    calc_value=.true., calc_derivative=.true.)

mu = 1._dp / natom * integral(L, ne * Hn)
mu_Hn = sum(Hn)/size(Hn)
print *, mu, mu_Hn
print *, "Summary of energies [a.u.]:"
print "('    Ts   = ', f14.8)", Ts
print "('    Een  = ', f14.8)", Een
print "('    Eee  = ', f14.8)", Eee
print "('    Exc  = ', f14.8)", Exc
print *, "   ---------------------"
print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV

print *, "Etot:", Etot, abs(Etot - Etot_conv)
print *, "mu:", mu, abs(mu - mu_conv)
print *, "mu_Hn:", mu_Hn, abs(mu_Hn - mu_conv)
print *, "abs(mu-mu_Hn):", abs(mu - mu_Hn)
call assert(abs(Etot - Etot_conv) < 1e-10_dp)
call assert(abs(mu - mu_conv) < 1e-5_dp)
call assert(abs(mu - mu_Hn) < 5e-5_dp)

end program
