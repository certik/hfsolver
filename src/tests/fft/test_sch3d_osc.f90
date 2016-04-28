program test_sch3d_osc

! Tests 1D harmonic oscillator ground state minimization (imaginary time
! propagation).

use types, only: dp
use constants, only: i_
use ofdft, only: read_pseudo
use ofdft_fft, only: free_energy, radial_potential_fourier, &
    reciprocal_space_vectors, free_energy_min, real2fourier, integral, &
    fourier2real, real_space_vectors, vtk_save, integralG
use utils, only: loadtxt, stop_error, assert, linspace, strfmt, init_random
use splines, only: spline3pars, iixmin, poly3, spline3ders
use interp3d, only: trilinear
use md, only: positions_fcc, positions_bcc
use linalg, only: eigvals, inv, solve
use sorting, only: argsort
implicit none
integer :: Ng
real(dp), allocatable :: G(:,:,:,:), G2(:,:,:), ne(:,:,:)
real(dp), allocatable :: Xn(:,:,:,:), Vn(:,:,:), r(:,:,:)
complex(dp), allocatable, dimension(:) :: lam
complex(dp), allocatable :: H(:,:), psi(:,:,:)
integer, allocatable :: idx(:)
integer, parameter :: nelec = 1
real(dp) :: L, dt, t, psi_norm
real(dp) :: omega, lambda
integer :: i, ai, aj, ak, bi, bj, bk, ci, cj, ck, a_idx, b_idx !, a, b
complex(dp) :: lam0

Ng = 8

L = 10._dp

allocate(G(Ng,Ng,Ng,3), G2(Ng,Ng,Ng), psi(Ng,Ng,Ng))
allocate(Xn(Ng, Ng, Ng, 3), Vn(Ng, Ng, Ng), r(Ng, Ng, Ng))
allocate(H(Ng**3, Ng**3))
allocate(lam(Ng**3), idx(Ng**3), ne(Ng,Ng,Ng))

call real_space_vectors(L, Xn)
call reciprocal_space_vectors(L, G, G2)
!omega = 1._dp
r = sqrt(sum((Xn-L/2)**2, dim=4))
!Vn = omega**2 * r**2 / 2
lambda = 0.2_dp
Vn = -exp(-lambda*r**2)

!a = 2
!b = 12
!Vn = -1/sqrt(a+(Xn-L/2)**b)

psi = 1 / L**3

dt = 1e-2_dp
print *, "dt =", dt

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

!    psi = psi * exp(-Vn*dt/2)
!    call real2fourier(psi, psiG)
!    psiG = psiG * exp(-G2*dt/2)
!    call fourier2real(psiG, psi)
!    psi = psi * exp(-Vn*dt/2)

!    ne = real(psi*conjg(psi), dp)
!    psi_norm = integral(L, ne)
!    psi = sqrt(nelec / psi_norm) * psi
!    ne = real(psi*conjg(psi), dp)
!    psi_norm = integral(L, ne)
!    print *, "norm of psi:", psi_norm

!    call real2fourier(psi, psiG)
!    E_tot = 1._dp/2 * integralG(G2*abs(psiG)**2, L) + integral(L, Vn*ne)
!    print *, "E_tot       =", E_tot
    print *, "E_tot_exact =", omega/2

end do

end program
