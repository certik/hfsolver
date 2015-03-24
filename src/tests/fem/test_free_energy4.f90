program test_free_energy4

! nuclear charge: 1 Gaussian
! electronic charge: 4 Gaussians
! calculation: single free energy evaluation

! This test uses FE and produces the same result as test_free_energy_fft2

use types, only: dp
use ofdft_fe, only: free_energy2, fe_data, initialize_fe
use ofdft_fft, only: reciprocal_space_vectors, radial_potential_fourier, &
    real2fourier
use constants, only: Ha2eV, pi, i_
use utils, only: loadtxt, assert, linspace
use splines, only: spline3pars, iixmin, poly3
use interp3d, only: trilinear
use feutils, only: quad_lobatto
use md, only: positions_fcc
use converged_energies, only: four_gaussians
use poisson3d_assembly, only: func2quad
use fe_mesh, only: quad2fe_3d, fe_eval_xyz
implicit none
real(dp) :: Eee, Een, Ts, Exc, Etot, Etot_conv
integer :: p, DOF, Nq
real(dp) :: Rcut, L, T_eV, T_au
integer, parameter :: natom = 4
real(dp) :: X(3, natom)
integer :: Nex, Ney, Nez
! This is for forces only:
real(dp) :: fen(3, natom)
real(dp), allocatable :: ne(:, :, :), R(:), Ven0G(:, :, :), fac(:, :, :), &
    G(:, :, :, :), G2(:, :, :), fullsol(:)
complex(dp), allocatable :: neG(:, :, :)
real(dp), allocatable, dimension(:, :, :, :) :: nq_pos
type(fe_data) :: fed
integer :: i, j, k, Ng


Rcut = 0.3_dp
p = 8
Nex = 8
Ney = 8
Nez = 8
L = 2
T_eV = 0.0862_dp
T_au = T_ev / Ha2eV
Nq = 9
call positions_fcc(X, L)
call free_energy2(real(natom, dp), L, Nex, Ney, Nez, p, T_au, nen, fne, &
        Nq, quad_lobatto, &
        Eee, Een, Ts, Exc, DOF)
Etot = Ts + Een + Eee + Exc
Etot_conv = sum(four_gaussians)
print *, "p =", p
print *, "DOF =", DOF
print *, "Rcut =", Rcut
print *, "T_au =", T_au
print *, "Summary of energies [a.u.]:"
print "('    Ts   = ', f14.8)", Ts
print "('    Een  = ', f14.8)", Een
print "('    Eee  = ', f14.8)", Eee
print "('    Exc  = ', f14.8)", Exc
print *, "   ---------------------"
print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV

print *, "Errors:"
print *, abs(Ts - four_gaussians(1))
print *, abs(Een - four_gaussians(2))
print *, abs(Eee - four_gaussians(3))
print *, abs(Exc - four_gaussians(4))
print *, abs(Etot - Etot_conv)
call assert(abs(Ts - four_gaussians(1)) < 1e-8_dp)
call assert(abs(Een - four_gaussians(2)) < 1e-8_dp)
call assert(abs(Eee - four_gaussians(3)) < 1e-8_dp)
call assert(abs(Exc - four_gaussians(4)) < 1e-8_dp)
call assert(abs(Etot - Etot_conv) < 1e-8_dp)

! ----------------------------------------------------------------------
! Forces calculation:

Ng = 80
allocate(Ven0G(Ng, Ng, Ng), ne(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng), fac(Ng, Ng, Ng), neG(Ng, Ng, Ng))


call initialize_fe(L, Nex, Ney, Nez, p, Nq, quad_lobatto, fed)
allocate(nq_pos(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(fullsol(maxval(fed%in)))
nq_pos = func2quad(fed%nodes, fed%elems, fed%xiq, fne)
call quad2fe_3d(fed%Ne, fed%Nb, fed%p, fed%jac_det, fed%wtq3, &
        fed%Sp, fed%Sj, fed%Sx, fed%phi_v, fed%in, fed%ib, &
        nq_pos, fullsol)
do i = 1, Ng
do j = 1, Ng
do k = 1, Ng
    ne(i, j, k) = fe_eval_xyz(fed%xin, fed%nodes, fed%elems, fed%in, &
        fullsol, [L, L, L]/Ng * ([i, j, k]-1))
end do
end do
end do

call reciprocal_space_vectors(L, G, G2)
allocate(R(40000))
R = linspace(1._dp/40000, 0.9_dp, 40000)
call radial_potential_fourier(R, 1*erf(6*R)/R, L, 1._dp, Ven0G)
call real2fourier(ne, neG)
fen = 0
do i = 1, natom
    fac = L**3*Ven0G*aimag(neG*exp(-i_ * &
        (G(:,:,:,1)*X(1,i) + G(:,:,:,2)*X(2,i) + G(:,:,:,3)*X(3,i))))
    fen(1, i) = sum(G(:,:,:,1)*fac)
    fen(2, i) = sum(G(:,:,:,2)*fac)
    fen(3, i) = sum(G(:,:,:,3)*fac)
end do

print *, "forces FE:"
print *, fen(:, 1)
print *, fen(:, 2)
print *, fen(:, 3)
print *, fen(:, 4)

contains

real(dp) function nen(x_, y_, z_) result(n)
real(dp), intent(in) :: x_, y_, z_
real(dp), parameter :: alpha = 6, Z = 1
real(dp) :: r2
integer :: i, a, b, c
n = 0
do i = 1, natom
    do a = -1, 1
    do b = -1, 1
    do c = -1, 1
        r2 = sum(([x_, y_, z_]-X(:, i)+[a, b, c]*L)**2)
        n = n - Z*alpha**3/pi**(3._dp/2)*exp(-alpha**2*r2)
    end do
    end do
    end do
end do
end function

real(dp) function fne(x, y, z) result(n)
real(dp), intent(in) :: x, y, z
real(dp), parameter :: alpha = 5, Z_ = 4
real(dp) :: r
r = sqrt((x-L/2)**2+(y-L/2)**2+(z-L/2)**2)
n = Z_*alpha**3/pi**(3._dp/2)*exp(-alpha**2*R**2)
end function

end program
