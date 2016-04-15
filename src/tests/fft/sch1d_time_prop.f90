program sch1d_time_prop

! Implements time propagation for a 1D Schroedinger equation.

use types, only: dp
use constants, only: i_, pi
use ofdft, only: read_pseudo
use ofdft_fft, only: free_energy, radial_potential_fourier, &
    reciprocal_space_vectors, free_energy_min, real2fourier, integral, &
    fourier2real, real_space_vectors, vtk_save
use utils, only: loadtxt, stop_error, assert, linspace, strfmt
use splines, only: spline3pars, iixmin, poly3, spline3ders
use interp3d, only: trilinear
use md, only: positions_fcc, positions_bcc
implicit none
integer :: Ng
real(dp), allocatable :: G(:), G2(:)
real(dp), allocatable :: ne(:)
real(dp), allocatable :: Xn(:), Vn(:)
complex(dp), allocatable, dimension(:) :: psi, psi2, psi3, psiG
complex(dp), allocatable, dimension(:) :: dpsi
real(dp) :: L
integer :: i
integer, parameter :: nelec = 1
real(dp) :: dt, psi_norm
integer :: u
real(dp) :: t

Ng = 32

L = 1._dp

allocate(ne(Ng))
allocate(G(Ng), G2(Ng), psi(Ng))
allocate(psi2(Ng), psi3(Ng), psiG(Ng))
allocate(dpsi(Ng))
allocate(Xn(Ng), Vn(Ng))

call real_space_vectors(L, Xn)
call reciprocal_space_vectors(L, G, G2)
Vn = 30*(Xn-L/2)**2
psi = gauss_x(Xn, 0.1_dp, L/2, -10._dp)

dt = 1e-4_dp
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

psi2 = psi

call real2fourier(psi2, psiG)
call fourier2real(-G2*psiG, dpsi)
psi = psi2 - i_*dt*(dpsi/2+Vn*psi2)

ne = real(psi*conjg(psi), dp)
psi_norm = integral(L, ne)
print *, "Initial norm of psi:", psi_norm
psi = sqrt(nelec / psi_norm) * psi
ne = real(psi*conjg(psi), dp)
psi_norm = integral(L, ne)
print *, "norm of psi:", psi_norm


do i = 1, 100
    t = t + dt
    print *, "iter =", i, "time =", t

    psi3 = psi2; psi2 = psi

    call real2fourier(psi2, psiG)
    call fourier2real(-G2*psiG, dpsi)

    psi = psi3 - 2*i_*dt*(dpsi/2+Vn*psi2)
    ne = real(psi*conjg(psi), dp)

    psi_norm = integral(L, ne)
    print *, "norm of psi:", psi_norm

    open(newunit=u, file="sch1d.txt", position="append", status="old")
    write(u, *) i, t, psi_norm
    close(u)

    open(newunit=u, file="sch1d_psi.txt", position="append", status="old")
    write(u, *) real(psi, dp), aimag(psi)
    close(u)

end do
print *, "Done"

contains

    pure function gauss_x(x, a, x0, k0) result(r)
    real(dp), intent(in) :: x(:), a, x0, k0
    complex(dp) :: r(size(x))
    r = (a * sqrt(pi))**(-0.5_dp) &
            * exp(-0.5_dp * ((x - x0) / a)**2 + i_ * x * k0)
    end function

end program
