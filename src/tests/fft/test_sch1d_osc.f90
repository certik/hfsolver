program test_sch1d_osc

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
real(dp), allocatable :: Xn(:), Vn(:)
complex(dp), allocatable, dimension(:) :: psi, psiG
real(dp) :: L
integer :: i
integer, parameter :: nelec = 1
real(dp) :: dt, psi_norm, E_tot, omega
integer :: u
real(dp) :: t

Ng = 16

L = 10._dp

allocate(ne(Ng))
allocate(G(Ng), G2(Ng), psi(Ng))
allocate(psiG(Ng))
allocate(Xn(Ng), Vn(Ng))

call real_space_vectors(L, Xn)
call reciprocal_space_vectors(L, G, G2)
omega = 1._dp
Vn = omega**2 * (Xn-L/2)**2 / 2
!Vn = reg_coulomb_x(Xn, 20._dp, L/2)
!psi = gauss_x(Xn, 0.1_dp, L/2, 400._dp)
psi = 1 / L

dt = 1e-2_dp
print *, "dt =", dt

open(newunit=u, file="sch1d.txt", status="replace")
close(u)

open(newunit=u, file="sch1d_grid.txt", status="replace")
write(u, *) Xn
write(u, *) Vn
close(u)

open(newunit=u, file="sch1d_psi.txt", status="replace")
write(u, *) real(psi, dp), aimag(psi)
close(u)


! Do first step by hand:
print *, "First step"
ne = real(psi*conjg(psi), dp)
psi_norm = integral(L, ne)
print *, "Initial norm of psi:", psi_norm
psi = sqrt(nelec / psi_norm) * psi
ne = real(psi*conjg(psi), dp)
psi_norm = integral(L, ne)
print *, "norm of psi:", psi_norm


t = 0

do i = 1, 600
    t = t + dt
    print *, "iter =", i, "time =", t

    psi = psi * exp(-Vn*dt/2)
    call real2fourier(psi, psiG)
    psiG = psiG * exp(-G2*dt/2)
    call fourier2real(psiG, psi)
    psi = psi * exp(-Vn*dt/2)

    ne = real(psi*conjg(psi), dp)
    psi_norm = integral(L, ne)
    psi = sqrt(nelec / psi_norm) * psi
    ne = real(psi*conjg(psi), dp)
    psi_norm = integral(L, ne)
    print *, "norm of psi:", psi_norm

    call real2fourier(psi, psiG)
    E_tot = 1._dp/2 * integralG(G2*abs(psiG)**2, L) + integral(L, Vn*ne)
    print *, "E_tot       =", E_tot
    print *, "E_tot_exact =", omega/2

    open(newunit=u, file="sch1d.txt", position="append", status="old")
    write(u, *) i, t, psi_norm, E_tot
    close(u)

    open(newunit=u, file="sch1d_psi.txt", position="append", status="old")
    write(u, *) real(psi, dp), aimag(psi)
    close(u)

end do
print *, "Done"

print *, E_tot - omega/2
call assert(abs(E_tot - omega/2) < 1e-9_dp)

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
