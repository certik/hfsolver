program test_free_energy2
use types, only: dp
use ofdft_fe, only: free_energy2
use ofdft, only: read_pseudo
use constants, only: Ha2eV, pi
use utils, only: loadtxt, stop_error, assert
use splines, only: spline3pars, iixmin, poly3, spline3ders
use interp3d, only: trilinear
use feutils, only: quad_gauss
use converged_energies, only: one_gaussian
implicit none
real(dp) :: Eh, Een, Ts, Exc, Etot
integer :: p, DOF, Nq
real(dp) :: Z, Ediff
real(dp), allocatable :: R(:), V(:), c(:, :)
real(dp), allocatable :: tmp(:), Vd(:), Vdd(:), density_en(:)
real(dp) :: Rcut, L, T_eV, T_au

call read_pseudo("H.pseudo.gaussian2", R, V, Z, Ediff)
allocate(tmp(size(R)), Vd(size(R)), Vdd(size(R)), density_en(size(R)))
call spline3ders(R, V, R, tmp, Vd, Vdd)
density_en = -(Vdd+2*Vd/R)/(4*pi)
allocate(c(0:4, size(R, 1)-1))
call spline3pars(R, density_en, [2, 2], [0._dp, 0._dp], c)
Rcut = R(size(R))
p = 5
L = 2
T_eV = 0.0862_dp
T_au = T_ev / Ha2eV
Nq = 20
call free_energy2(1._dp, L, 4, 4, 4, p, T_au, nen_splines, ne, &
        Nq, quad_gauss, &
        Eh, Een, Ts, Exc, DOF)
Etot = Eh + Een + Ts + Exc
print *, "Summary of energies [a.u.]:"
print "('    Ts   = ', f14.8)", Ts
print "('    Een  = ', f14.8)", Een
print "('    Eee  = ', f14.8)", Eh
print "('    Exc  = ', f14.8)", Exc
print *, "   ---------------------"
print "('    Etot = ', f14.8, ' a.u.')", Etot
call assert(abs(Ts - one_gaussian(1)) < 1e-8_dp)
call assert(abs(Een - (-3.77207363_dp)) < 3e-3_dp)
call assert(abs(Eh - one_gaussian(3)) < 3e-4_dp)
call assert(abs(Exc - one_gaussian(4)) < 1e-8_dp)

contains

real(dp) function nen_splines(x_, y_, z_) result(n)
real(dp), intent(in) :: x_, y_, z_
real(dp) :: r_
integer :: ip
! One atom in the center:
r_ = sqrt((x_+L/64-L/2)**2+(y_-L/2)**2+(z_-L/2)**2)
if (r_ >= Rcut) then
    n = 0
    return
end if
if (r_ <= 1e-4_dp) r_ = 1e-4_dp
ip = 0
ip = iixmin(r_, R, ip)
n = poly3(r_, c(:, ip))
! FIXME: We need to be using "rho", and flip the sign here:
n = -n
end function

real(dp) function ne(x, y, z) result(n)
real(dp), intent(in) :: x, y, z
real(dp), parameter :: alpha = 5, Z_ = 1
real(dp) :: r
r = sqrt((x-L/2)**2+(y-L/2)**2+(z-L/2)**2)
n = Z_*alpha**3/pi**(3._dp/2)*exp(-alpha**2*R**2)
!n = 1
end function

end program
