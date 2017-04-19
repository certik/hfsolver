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
integer :: u, u2
real(dp) :: t

Ng = 1024

L = 8

allocate(ne(Ng))
allocate(G(Ng), G2(Ng), psi(Ng))
allocate(psiG(Ng))
allocate(Xn(Ng), Vn(Ng))

call real_space_vectors(L, Xn)
call reciprocal_space_vectors(L, G, G2)

open(newunit=u, file="sch1d_grid.txt", status="replace")
write(u, *) Xn
Vn = gaussian_potential(Xn, 2._dp, L/2)
write(u, *) Vn
psi = gaussian_density(Xn, 2._dp, L/2)
write(u, *) real(psi, dp)

! Solve Poisson
call real2fourier(psi, psiG)
psiG(1) = 0; psiG(2:) = 4*pi*psiG(2:) / G2(2:)
call fourier2real(psiG, Vn)

write(u, *) Vn
close(u)

psi = 1 / L

dt = 1e-2_dp
print *, "dt =", dt

open(newunit=u, file="sch1d.txt", status="replace")

open(newunit=u2, file="sch1d_psi.txt", status="replace")
write(u2, *) real(psi, dp), aimag(psi)


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

    write(u, *) i, t, psi_norm, E_tot
    write(u2, *) real(psi, dp), aimag(psi)

end do
print *, "Done"

close(u)
close(u2)

print *, E_tot - omega/2
call assert(abs(E_tot - omega/2) < 1e-9_dp)

contains

    pure function gaussian_density(x, alpha, x0) result(n)
    ! Returns a Gaussian charge.
    !
    ! This function returns the particle density `n` corresponding to the
    ! following radial potential `V`:
    !     V = -erf(alpha*r)/r
    ! by solving the radial Poisson equation (r*V(r))'' = -4*pi*r*n(r), e.g. in
    ! SymPy:
    !
    !     >>> var("alpha r", positive=True)
    !     >>> V = -erf(alpha*r)/r
    !     >>> n = (r*V).diff(r, 2) / (-4*pi*r)
    !     >>> n
    !     -alpha**3*exp(-alpha**2*r**2)/pi**(3/2)
    !
    ! Note that in 1D the potential `V` as calculated from the 1D Poisson
    ! equation V''(x) = -4*pi*n(x) is equal to:
    !
    !     >>> simplify(integrate(integrate(-4*pi*n, r), r))
    !     2*alpha**2*r*erf(alpha*r) + 2*alpha*exp(-alpha**2*r**2)/sqrt(pi)
    !
    ! The last potential is returned by gaussian_potential().
    real(dp), intent(in) :: x(:), alpha, x0
    real(dp) :: n(size(x))
    n = -alpha**3 * exp(-alpha**2*(x-x0)**2)/pi**(3._dp/2)
    end function

    pure function gaussian_potential(x, alpha, x0) result(V)
    ! Returns a Gaussian potential. See gaussian_density() for details.
    real(dp), intent(in) :: x(:), alpha, x0
    real(dp) :: V(size(x)), r(size(x))
    r = abs(x-x0)
    V = 2*alpha**2*r*erf(alpha*r) + 2*alpha*exp(-alpha**2*r**2)/sqrt(pi)
    end function

end program
