program test_sch1d_osc_eigen

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
use linalg, only: eigvals
implicit none
integer :: Ng
real(dp), allocatable :: G(:), G2(:)
real(dp), allocatable :: ne(:)
real(dp), allocatable :: Xn(:), Vn(:)
complex(dp), allocatable, dimension(:) :: psi, psiG, lam
complex(dp), allocatable :: H(:,:)
real(dp) :: L
real(dp) :: dt, psi_norm, E_tot, omega
integer :: i, j, k

Ng = 16

L = 10._dp

allocate(ne(Ng))
allocate(G(Ng), G2(Ng), psi(Ng))
allocate(psiG(Ng))
allocate(Xn(Ng), Vn(Ng))
allocate(H(Ng, Ng))
allocate(lam(Ng))

call real_space_vectors(L, Xn)
call reciprocal_space_vectors(L, G, G2)
omega = 1._dp
Vn = omega**2 * (Xn-L/2)**2 / 2

call real2fourier(Vn, psiG)

do j = 1, Ng
do i = 1, Ng
    k = i-j+1
    if (k < 1) k = k + Ng
    H(i,j) = psiG(k)
end do
end do

do i = 1, Ng
    H(i,i) = H(i,i) + G2(i)/2
end do

lam = eigvals(H)
print *, "Eigenvalues:"
do i = 1, Ng
    print *, i, lam(i)
end do

psi = 1 / L

dt = 0.01_dp
psi = psi * exp(-Vn*dt/2)
call real2fourier(psi, psiG)
psiG = psiG * exp(-G2*dt/2)
call fourier2real(psiG, psi)
psi = psi * exp(-Vn*dt/2)

ne = real(psi*conjg(psi), dp)
psi_norm = integral(L, ne)
print *, "norm of psi:", psi_norm

call real2fourier(psi, psiG)
E_tot = 1._dp/2 * integralG(G2*abs(psiG)**2, L) + integral(L, Vn*ne)
print *, "E_tot       =", E_tot
print *, "E_tot_exact =", omega/2

print *, E_tot - omega/2
call assert(abs(E_tot - omega/2) < 1e-9_dp)

end program
