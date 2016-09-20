program test_functional_derivative_fft

! The same as test_functional_derivative but using FFT

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
use md, only: positions_fcc
use converged_energies, only: four_gaussians_min
implicit none
real(dp) :: Eee, Een, Ts, Exc, Etot, Etot_conv
integer :: Ng
real(dp) :: Z
real(dp), allocatable :: R(:), G(:, :, :, :), G2(:, :, :)
real(dp), allocatable :: ne(:, :, :), Hn(:, :, :)
real(dp), allocatable :: Ven0G(:, :, :), current(:,:,:,:), Xn(:,:,:,:)
real(dp) :: V0
complex(dp), allocatable, dimension(:,:,:) :: VenG, psi, psi2, psi3, psiG, tmp
complex(dp), allocatable, dimension(:,:,:,:) :: dpsi
real(dp) :: L, T_eV, T_au
integer :: i, j
integer, parameter :: natom = 4
real(dp) :: X(3, natom), alpha_nen, mu, dt, psi_norm
integer :: cg_iter
real(dp) :: E0, t, omega, current_avg(3), conductivity

Ng = 32

L = 2
T_eV = 0.0862_dp
T_au = T_ev / Ha2eV

alpha_nen = 6

Z = 1

allocate(Ven0G(Ng, Ng, Ng), VenG(Ng, Ng, Ng), ne(Ng, Ng, Ng), Hn(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng), psi(Ng, Ng, Ng))
allocate(R(40000), psi2(Ng, Ng, Ng), psi3(Ng, Ng, Ng), psiG(Ng, Ng, Ng))
allocate(dpsi(Ng, Ng, Ng, 3), current(Ng, Ng, Ng, 3), tmp(Ng, Ng, Ng))
allocate(Xn(Ng, Ng, Ng, 3))
R = linspace(1._dp/40000, 0.9_dp, 40000)
print *, "Radial nuclear potential FFT"
call radial_potential_fourier(R, Z*erf(alpha_nen*R)/R, L, Z, Ven0G, V0)
print *, "    Done."

call real_space_vectors([L, L, L], Xn)
call reciprocal_space_vectors([L, L, L], G, G2)
call positions_fcc(X, L)
VenG = 0
do i = 1, natom
    VenG = VenG - Ven0G * exp(-i_ * &
        (G(:,:,:,1)*X(1,i) + G(:,:,:,2)*X(2,i) + G(:,:,:,3)*X(3,i)))
end do

! Minimization

ne = natom / L**3
call free_energy_min(real(natom, dp), natom, L, G2, T_au, VenG, ne, &
    1e-12_dp, Eee, Een, Ts, Exc, Etot, cg_iter)

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
call assert(abs(Ts - four_gaussians_min(1)) < 1e-7_dp)
call assert(abs(Een - four_gaussians_min(2)) < 1e-7_dp)
call assert(abs(Eee - four_gaussians_min(3)) < 1e-8_dp)
call assert(abs(Exc - four_gaussians_min(4)) < 1e-8_dp)
call assert(abs(Etot - Etot_conv) < 1e-10_dp)

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
print *, "dt =", dt

! Do first step by hand:
print *, "First step"
psi = sqrt(ne)

t = 0

psi2 = psi
psi = psi2 - i_*dt*Hn*psi2

ne = real(psi*conjg(psi), dp)
psi_norm = integral(L, ne)
print *, "Initial norm of psi:", psi_norm
psi = sqrt(natom / psi_norm) * psi
ne = real(psi*conjg(psi), dp)
psi_norm = integral(L, ne)
print *, "norm of psi:", psi_norm

E0 = 1e-3_dp
omega = 0.05

do i = 1, 10
    t = t + dt
    print *, "iter =", i, "time =", t
    psi3 = psi2; psi2 = psi
    psi = psi3 - 2*i_*dt*(Hn + E0*Xn(:,:,:,1)*sin(omega*t))*psi2
    ne = real(psi*conjg(psi), dp)
    call real2fourier(psi, psiG)
    psiG(1,1,1) = 0
    do j = 1, 3
        call fourier2real(i_*G(:,:,:,j)*psiG, dpsi(:,:,:,j))
        tmp = (conjg(psi)*dpsi(:,:,:,j)-psi*conjg(dpsi(:,:,:,j))) / (2*natom*i_)
        if (maxval(abs(aimag(tmp))) > 1e-12_dp) then
            print *, "INFO: current  max imaginary part:", maxval(aimag(tmp))
        end if
        current(:,:,:,j) = real(tmp, dp)
    end do

    psi_norm = integral(L, ne)
    print *, "norm of psi:", psi_norm

    call free_energy(L, G2, T_au, VenG, ne, Eee, Een, Ts, Exc, Etot, Hn, &
        calc_value=.true., calc_derivative=.true.)
    Etot = Ts + Een + Eee + Exc
    print *, "Summary of energies [a.u.]:"
    print "('    Ts   = ', f14.8)", Ts
    print "('    Een  = ', f14.8)", Een
    print "('    Eee  = ', f14.8)", Eee
    print "('    Exc  = ', f14.8)", Exc
    print *, "   ---------------------"
    print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV


    do j = 1, 3
        current_avg(j) = integral(L, current(:, :, :, j))/L**3
    end do
    print *, "E field along the 'x' direction =", E0*sin(omega*t)
    print *, "average current =", current_avg
    print *, "current normalized =", current_avg / current_avg(1)
    conductivity = current_avg(1) / E0*sin(omega*t)
    print *, "conductivity along the 'x' direction =", conductivity

end do
print *, "Done"

end program
