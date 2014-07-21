program test_minimization
use types, only: dp
use ofdft_fe, only: free_energy_min
use ofdft, only: read_pseudo
use constants, only: Ha2eV, pi
use utils, only: loadtxt, assert
use splines, only: spline3pars, iixmin, poly3
use interp3d, only: trilinear
use feutils, only: quad_gauss
implicit none
real(dp) :: Eh, Een, Ts, Exc, Etot
integer :: p, DOF !, Nx, Ny, Nz, u
real(dp) :: Z, Ediff
real(dp), allocatable :: R(:), V(:), D(:, :), c(:, :), values(:, :, :)
real(dp) :: Rcut, L, T_eV, T_au

!open(newunit=u, file="plots/Ven_reg128.txt", status="old")
!read(u, *) Nx, Ny, Nz
!allocate(values(Nx, Ny, Nz))
!read(u, *) values
!close(u)
!call read_pseudo("H.pseudo", R, V, Z, Ediff)
!V = -V
allocate(R(1), V(1)); Z=1; Ediff=0
!call loadtxt("Venr.txt", D)
!allocate(c(0:4, size(D, 1)-1))
!call spline3pars(D(:, 1), D(:, 2), [2, 2], [0._dp, 0._dp], c)
Rcut = R(size(R))
Rcut = 0.3_dp
p = 2
L = 2
T_eV = 0.0862_dp
T_au = T_ev / Ha2eV
call free_energy_min(1._dp, L, 2, 2, 2, p, T_au, nen, ne, &
    5, quad_gauss, 1e-9_dp, &
    Eh, Een, Ts, Exc, DOF)
Etot = Ts + Een + Eh + Exc
print *, "p =", p
print *, "DOF =", DOF
print *, "Rcut =", Rcut
print *, "T_au =", T_au
print *, "Summary of energies [a.u.]:"
print "('    Ts   = ', f14.8)", Ts
print "('    Een  = ', f14.8)", Een
print "('    Eee  = ', f14.8)", Eh
print "('    Exc  = ', f14.8)", Exc
print *, "   ---------------------"
print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV

! These values are specific to the given Nq, p, Nx, Ny, Nz above, and
! self-consistency 1e-9. The energy components (Ts, Een, Eh, Exc) have lower
! accuracy than the total energy.
call assert(abs(Ts-1.5262_dp) < 1e-4_dp)
call assert(abs(Een-(-1.4563_dp)) < 1e-4_dp)
call assert(abs(Eh-0.0981_dp) < 1e-4_dp)
call assert(abs(Exc-(-0.5481_dp)) < 1e-4_dp)
call assert(abs(Etot-(-0.38016292_dp)) < 1e-8_dp)

contains

real(dp) function nen(x, y, z_) result(n)
real(dp), intent(in) :: x, y, z_
real(dp), parameter :: alpha = 12
real(dp) :: r
r = sqrt((x-L/2)**2+(y-L/2)**2+(z_-L/2)**2)
! This density:
n = -Z*alpha**3/pi**(3._dp/2)*exp(-alpha**2*R**2)
! Corresponds to the potential:
!V = -Z*erf(alpha*R)/R
end function

real(dp) function Ven_splines(x_, y_, z_) result(V)
real(dp), intent(in) :: x_, y_, z_
real(dp) :: r
integer :: ip
! One atom in the center:
r = sqrt(x_**2+y_**2+z_**2)
if (r >= 1) r = 1
ip = 0
ip = iixmin(r, D(:, 1), ip)
V = poly3(r, c(:, ip))
end function

real(dp) function Ven_interp3d(x, y, z) result(V)
real(dp), intent(in) :: x, y, z
V = trilinear([x, y, z], [-L/2, -L/2, -L/2], [L/2, L/2, L/2], &
        values)
end function

real(dp) function ne(x, y, z) result(n)
real(dp), intent(in) :: x, y, z
real(dp), parameter :: alpha = 1, Z_ = 1
real(dp) :: r
r = sqrt((x-L/2)**2+(y-L/2)**2+(z-L/2)**2)
n = Z_*alpha**3/pi**(3._dp/2)*exp(-alpha**2*R**2)
n = 1
end function

end program
