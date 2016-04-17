program test_sch1d_osc2

! Tests 1D harmonic oscillator ground state minimization (imaginary time
! propagation).

use types, only: dp
use constants, only: i_, pi
use ofdft, only: read_pseudo
use ofdft_fft, only: free_energy, radial_potential_fourier, &
    reciprocal_space_vectors, free_energy_min, real2fourier, integral, &
    fourier2real, real_space_vectors, vtk_save, integralG
use utils, only: loadtxt, stop_error, assert, linspace, strfmt
use splines, only: spline3pars, iixmin, poly3, spline3ders
use interp3d, only: trilinear
use md, only: positions_fcc, positions_bcc
implicit none
integer :: Ng
real(dp), allocatable :: G(:), G2(:)
real(dp), allocatable :: ne(:)
real(dp), allocatable :: Xn(:), Vn(:), current(:)
complex(dp), allocatable, dimension(:) :: psi, psiG, tmp, dpsi
real(dp) :: L
integer :: i
integer, parameter :: nelec = 1
real(dp) :: dt, psi_norm, E_tot, current_avg, Ex, E0, td, tw, a, V0, b
integer :: u, u2, Nmult
real(dp) :: t

Nmult = 16

Ng = 64 * Nmult

L = 10._dp * Nmult

allocate(ne(Ng))
allocate(G(Ng), G2(Ng), psi(Ng))
allocate(psiG(Ng), tmp(Ng), dpsi(Ng), current(Ng))
allocate(Xn(Ng), Vn(Ng))

call real_space_vectors(L, Xn)
call reciprocal_space_vectors(L, G, G2)
a = 2
b = 12
Vn = -1/sqrt(a+(Xn-L/2)**b)

open(newunit=u2, file="sch1d_psi2.txt", status="replace")
write(u2, *) real(psi, dp), aimag(psi)

psi = 1 / L


open(newunit=u, file="sch1d_grid2.txt", status="replace")
write(u, *) Xn
write(u, *) Vn
close(u)


open(newunit=u, file="sch1d2.txt", status="replace")

open(newunit=u2, file="sch1d_psi2.txt", status="replace")
write(u2, *) real(psi, dp), aimag(psi)

dt = 1e-3_dp
E0 = 0.001_dp
td = 0.2_dp
tw = 0.04_dp

print *, "dt =", dt

t = 0
do i = 1, 100
    t = t + dt
    print *, "iter =", i, "time =", t
    Ex = E0 * exp(-(t-td)**2/(2*tw**2)) / (sqrt(2*pi)*tw)
!    Ex = 0

    psi = psi * exp(-i_*(Vn+Xn*Ex)*dt/2)
    call real2fourier(psi, psiG)
    psiG = psiG * exp(-i_*G2*dt/2)
    call fourier2real(psiG, psi)
    psi = psi * exp(-i_*(Vn+Xn*Ex)*dt/2)

    ne = real(psi*conjg(psi), dp)
    psi_norm = integral(L, ne)
    print *, "norm of psi:", psi_norm

    call real2fourier(psi, psiG)
    E_tot = 1._dp/2 * integralG(G2*abs(psiG)**2, L) &
        + integral(L, (Vn+Xn*Ex)*ne)
    print *, "E_tot       =", E_tot

    call fourier2real(i_*G*psiG, dpsi)
    tmp = (conjg(psi)*dpsi-psi*conjg(dpsi)) / (2*i_)
    current = real(tmp, dp)
    call assert(maxval(abs(aimag(tmp))) < epsilon(1._dp))
    current_avg = integral(L, current)/L
    print *, "average current =", current_avg

    write(u, *) i, t, psi_norm, E_tot, current_avg, Ex
    write(u2, *) real(psi, dp), aimag(psi)

end do

close(u)
close(u2)

contains

    pure function gauss_x(x, a, x0, k0) result(r)
    real(dp), intent(in) :: x(:), a, x0, k0
    complex(dp) :: r(size(x))
    r = (a * sqrt(pi))**(-0.5_dp) &
            * exp(-0.5_dp * ((x - x0) / a)**2 + i_ * x * k0)
    end function

    pure function reg_coulomb_x(x, alpha, x0) result(V)
    real(dp), intent(in) :: x(:), alpha, x0
    real(dp) :: V(size(x)), r(size(x))
    r = abs(x-x0)
    where (r < 1e-6_dp)
        V = -2*alpha/sqrt(pi)
    else where
        V = -erf(alpha*r)/r
    end where
    end function

end program
