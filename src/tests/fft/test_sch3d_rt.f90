program test_sch3d_rt

! 3D RT

use types, only: dp
use constants, only: i_, pi
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
real(dp), allocatable :: Xn(:,:,:,:), Vn(:,:,:), r(:,:,:), current(:,:,:,:)
complex(dp), allocatable :: psi(:,:,:), psiG(:,:,:), dpsi(:,:,:,:), tmp(:,:,:)
integer, parameter :: nelec = 1
real(dp) :: L, dt, t, psi_norm, E_tot, current_avg(3)
real(dp) :: lambda, E0, Ex, td, tw
integer :: i, j, u, a, b, Nmult

Nmult = 8

Ng = 8 * Nmult

L = 10._dp * Nmult

allocate(G(Ng,Ng,Ng,3), G2(Ng,Ng,Ng), psi(Ng,Ng,Ng))
allocate(Xn(Ng, Ng, Ng, 3), Vn(Ng, Ng, Ng), r(Ng, Ng, Ng))
allocate(psiG(Ng,Ng,Ng))
allocate(ne(Ng,Ng,Ng), tmp(Ng,Ng,Ng))
allocate(dpsi(Ng,Ng,Ng,3), current(Ng,Ng,Ng,3))

call real_space_vectors(L, Xn)
call reciprocal_space_vectors(L, G, G2)
!lambda = 0.2_dp
!Vn = -exp(-lambda*r**2)
a = 2
b = 12
r = sqrt(sum((Xn-L/2)**2, dim=4))
Vn = -3/sqrt(a+r**b)

ne = 1 / L**3
psi = sqrt(ne)

dt = 0.1_dp
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

do i = 1, 100
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

end do

print *
print *, "Real Time Propagation:"

open(newunit=u, file="cond.txt", status="replace")

dt = 1e-2_dp
E0 = 0.001_dp
td = 0.2_dp
tw = 0.04_dp
t = 0

do i = 1, 5000
    t = t + dt
    print *, "iter =", i, "time =", t
    Ex = E0 * exp(-(t-td)**2/(2*tw**2)) / (sqrt(2*pi)*tw)

    psi = psi * exp(-i_*(Vn+Xn(:,:,:,1)*Ex)*dt/2)
    call real2fourier(psi, psiG)
    psiG = psiG * exp(-i_*G2*dt/2)
    call fourier2real(psiG, psi)
    psi = psi * exp(-i_*(Vn+Xn(:,:,:,1)*Ex)*dt/2)

    ne = real(psi*conjg(psi), dp)
    psi_norm = integral(L, ne)
    print *, "norm of psi:", psi_norm

    ! Total energy
    call real2fourier(psi, psiG)
    E_tot = 1._dp/2 * integralG(G2*abs(psiG)**2, L) + integral(L, Vn*ne)
    print *, "E_tot       =", E_tot

    ! Current
    do j = 1, 3
        call fourier2real(i_*G(:,:,:,j)*psiG, dpsi(:,:,:,j))
        tmp = (conjg(psi)*dpsi(:,:,:,j)-psi*conjg(dpsi(:,:,:,j))) / (2*nelec*i_)
        call assert(maxval(abs(aimag(tmp))) < epsilon(1._dp))
        current(:,:,:,j) = real(tmp, dp)
    end do
    do j = 1, 3
        current_avg(j) = integral(L, current(:, :, :, j))/L**3
    end do
    print *, "average current X =", current_avg(1)
    print *, "average current Y =", current_avg(2)
    print *, "average current Z =", current_avg(3)

    write(u, *) i, t, Ex, E_tot, psi_norm, current_avg

end do

close(u)

end program
