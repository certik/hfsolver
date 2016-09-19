program ofmd_rt

! Real time (RT) propagation using OFMD.

use types, only: dp
use constants, only: i_, pi
use ofdft, only: read_pseudo
use ofdft_fft, only: free_energy, radial_potential_fourier, &
    reciprocal_space_vectors, free_energy_min, real2fourier, integral, &
    fourier2real, real_space_vectors
use constants, only: Ha2eV
use utils, only: loadtxt, stop_error, assert, linspace
use splines, only: spline3pars, iixmin, poly3, spline3ders
use interp3d, only: trilinear
use md, only: positions_fcc
implicit none
real(dp) :: Eee, Een, Ts, Exc, Etot
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
integer, parameter :: natom = 1
real(dp) :: X(3, natom), alpha_nen, mu, dt, psi_norm
integer :: cg_iter
real(dp) :: E0, t, current_avg(3), Ex, td, tw

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

call real_space_vectors(L, Xn)
call reciprocal_space_vectors(L, G, G2)
!call positions_fcc(X, L)
! Put the single atom into the middle of the box [0, L]^3.
X(:, 1) = [L/2, L/2, L/2]
VenG = 0
do i = 1, natom
    VenG = VenG - Ven0G * exp(-i_ * &
        (G(:,:,:,1)*X(1,i) + G(:,:,:,2)*X(2,i) + G(:,:,:,3)*X(3,i)))
end do

! Minimization

ne = natom / L**3
call free_energy_min(real(natom, dp), natom, L, G2, T_au, VenG, ne, &
    1e-12_dp, Eee, Een, Ts, Exc, Etot, cg_iter)

print *, "Summary of energies [a.u.]:"
print "('    Ts   = ', f14.8)", Ts
print "('    Een  = ', f14.8)", Een
print "('    Eee  = ', f14.8)", Eee
print "('    Exc  = ', f14.8)", Exc
print *, "   ---------------------"
print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV

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
! This is an automatic way to get the time step, but for now we leave it
! commented out and set the time step by hand below:
!dt = 1/maxval(abs(Hn)) / 10 ! set dt 10x smaller than the limit

! Set the time step by hand:
dt = 1e-2_dp

E0 = 0.001_dp
td = 0.2_dp
tw = 0.04_dp
print *, "dt =", dt

! Do first step by hand:
print *, "First step"
psi = sqrt(ne)

t = 0

do i = 1, 10
    t = t + dt
    print *, "iter =", i, "time =", t
    Ex = E0 * exp(-(t-td)**2/(2*tw**2)) / (sqrt(2*pi)*tw)

    psi = psi * exp(-i_*(Hn+Xn(:,:,:,1)*Ex)*dt/2)
! There is no Laplace operator, so this is skipped:
!    call real2fourier(psi, psiG)
!    psiG = psiG * exp(-i_*G2*dt/2)
!    call fourier2real(psiG, psi)
    psi = psi * exp(-i_*(Hn+Xn(:,:,:,1)*Ex)*dt/2)

    ne = real(psi*conjg(psi), dp)
    call real2fourier(psi, psiG)
    ! This is probably not needed, so it's commented out:
    !psiG(1,1,1) = 0
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
    print *, "average current =", current_avg
    print *, "current normalized =", current_avg / current_avg(1)

end do
print *, "Done"

end program
