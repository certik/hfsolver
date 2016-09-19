program ofmd_time_propagation

! Implements time propagation for a given OFMD MD step.

use types, only: dp
use constants, only: i_, pi
use ofdft, only: read_pseudo
use ofdft_fft, only: free_energy, radial_potential_fourier, &
    reciprocal_space_vectors, free_energy_min, real2fourier, integral, &
    fourier2real, real_space_vectors, vtk_save
use constants, only: Ha2eV
use utils, only: loadtxt, stop_error, assert, linspace, strfmt
use splines, only: spline3pars, iixmin, poly3, spline3ders
use interp3d, only: trilinear
use md, only: positions_fcc, positions_bcc
implicit none
real(dp) :: Eee, Een, Ts, Exc, Etot, Etot_conv
integer :: Ng
real(dp) :: Z
real(dp), allocatable :: R(:), G(:, :, :, :), G2(:, :, :), Ven_rad(:)
real(dp), allocatable :: ne(:, :, :), Hn(:, :, :)
real(dp), allocatable :: Ven0G(:, :, :), current(:,:,:,:), Xn(:,:,:,:)
real(dp) :: V0
complex(dp), allocatable, dimension(:,:,:) :: VenG, psi, psi2, psi3, psiG, tmp
complex(dp), allocatable, dimension(:,:,:,:) :: dpsi
real(dp) :: L, T_eV, T_au
integer :: i, j
integer, parameter :: natom = 128
real(dp) :: mu, dt, psi_norm
real(dp), allocatable :: X(:, :)
integer :: cg_iter, u
real(dp) :: Ex, E0, t, current_avg(3), conductivity, td, tw, Ediff

Ng = 64

L = 8.1049178668765851_dp
T_eV = 34.5_dp
T_au = T_ev / Ha2eV

allocate(Ven0G(Ng, Ng, Ng), VenG(Ng, Ng, Ng), ne(Ng, Ng, Ng), Hn(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng), psi(Ng, Ng, Ng))
allocate(psi2(Ng, Ng, Ng), psi3(Ng, Ng, Ng), psiG(Ng, Ng, Ng))
allocate(dpsi(Ng, Ng, Ng, 3), current(Ng, Ng, Ng, 3), tmp(Ng, Ng, Ng))
allocate(Xn(Ng, Ng, Ng, 3), X(3, natom))
call read_pseudo("D.pseudo", R, Ven_rad, Z, Ediff)
print *, "Radial nuclear potential FFT"
call radial_potential_fourier(R, Ven_rad, L, Z, Ven0G, V0)
print *, "    Done."

call real_space_vectors([L, L, L], Xn)
call reciprocal_space_vectors([L, L, L], G, G2)
call positions_bcc(X, L)
!call loadtxt("../pos.txt", X)
call assert(size(X, 1) == 3)
call assert(size(X, 2) == natom)
call assert(all(X > 0))
call assert(all(X < L))
VenG = 0
do i = 1, natom
    VenG = VenG - Ven0G * exp(-i_ * &
        (G(:,:,:,1)*X(1,i) + G(:,:,:,2)*X(2,i) + G(:,:,:,3)*X(3,i)))
end do

! Minimization

ne = natom / L**3
call free_energy_min(real(natom, dp), natom, L, G2, T_au, VenG, ne, &
    2e-10_dp, Eee, Een, Ts, Exc, Etot, cg_iter)

Etot_conv = -204.88460758_dp
print *, "Summary of energies [a.u.]:"
print "('    Ts   = ', f14.8)", Ts
print "('    Een  = ', f14.8)", Een
print "('    Eee  = ', f14.8)", Eee
print "('    Exc  = ', f14.8)", Exc
print *, "   ---------------------"
print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV
print "('    Etot/atom = ', f14.8, ' a.u. = ', f14.8, ' eV')", &
    Etot/natom, Etot*Ha2eV / natom
print *, "Errors:"
print *, abs(Etot - Etot_conv)
!call assert(abs(Etot - Etot_conv) < 3e-3_dp)

call free_energy(L, G2, T_au, VenG, ne, Eee, Een, Ts, Exc, Etot, Hn, &
    calc_value=.true., calc_derivative=.true.)

mu = 1._dp / natom * integral(L, ne * Hn)
print *, "mu    = ", mu
mu = sum(Hn)/size(Hn)
print *, "mu_Hn = ", mu
print *, "max(abs(H-mu)) = ", maxval(abs(Hn - mu))
call assert(all(ne > 0))

print *
print *, "------------------------------------------------------------------"
print *, "Propagation"

! Propagate

print *, "E_max =", maxval(abs(Hn)), "; dt <", 1/maxval(abs(Hn))
dt = 1/maxval(abs(Hn)) / 20 ! set dt 10x smaller than the limit
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

E0 = 0.01_dp


open(newunit=u, file="cond.txt", status="replace")
close(u)

td = 0.2_dp
tw = 0.04_dp

do i = 1, 1000
    t = t + dt
    print *, "iter =", i, "time =", t
    print *, "dt     =", dt
    print *, "dt max =", 1/maxval(abs(Hn))
    if (dt > 1/maxval(abs(Hn))) &
        call stop_error("Time step is too large (dt > dt_max).")
    psi3 = psi2; psi2 = psi
    Ex = E0 * exp(-(t-td)**2/(2*tw**2)) / (sqrt(2*pi)*tw)
    psi = psi3 - 2*i_*dt*(Hn + Xn(:,:,:,1)*Ex)*psi2
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
    print "('    Etot/atom = ', f14.8, ' a.u. = ', f14.8, ' eV')", &
        Etot/natom, Etot*Ha2eV / natom


    do j = 1, 3
        current_avg(j) = integral(L, current(:, :, :, j))/L**3
    end do
    print *, "E field along the 'x' direction =", Ex
    print *, "average current =", current_avg
    print *, "current normalized =", current_avg / current_avg(1)
    conductivity = current_avg(1) / E0
    print *, "conductivity along the 'x' direction =", conductivity

    open(newunit=u, file="cond.txt", position="append", status="old")
    write(u, *) i, t, Etot*Ha2eV / natom, Ex, conductivity, current_avg, &
        psi_norm, dt, 1/maxval(abs(Hn))
    close(u)

    !call vtk_save("data/iter" // strfmt('(i0.6)', i) // ".vtk", Xn, ne)

end do
print *, "Done"

end program
