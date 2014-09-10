program test_free_energy3

! This test uses FE and produces the same result as test_free_energy_fft

use types, only: dp
use ofdft_fe, only: free_energy2
use constants, only: Ha2eV, pi
use utils, only: loadtxt, assert
use splines, only: spline3pars, iixmin, poly3
use interp3d, only: trilinear
use feutils, only: quad_gauss
implicit none
real(dp) :: Eee, Een, Ts, Exc, Etot
integer :: p, DOF, Nq
real(dp) :: Z
real(dp) :: Rcut, L, T_eV, T_au

Z = 1
Rcut = 0.3_dp
p = 8
L = 2
T_eV = 0.0862_dp
T_au = T_ev / Ha2eV
Nq = 25
call free_energy2(1._dp, L, 8, 8, 8, p, T_au, nen, ne, &
        Nq, quad_gauss, &
        Eee, Een, Ts, Exc, DOF)
Etot = Ts + Een + Eee + Exc
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
print *, abs(Ts - (10.61904507_dp))
print *, abs(Een - (-2.92172113_dp))
print *, abs(Eee - (1.30109500_dp))
print *, abs(Exc - (-1.43805890_dp))
print *, abs(Etot - (7.56036004_dp))
call assert(abs(Ts - (10.61904507_dp)) < 1e-8_dp)
call assert(abs(Een - (-2.92172113_dp)) < 1e-8_dp)
call assert(abs(Eee - (1.30109500_dp)) < 1e-8_dp)
call assert(abs(Exc - (-1.43805890_dp)) < 1e-8_dp)
call assert(abs(Etot - (7.56036004_dp)) < 1e-8_dp)

contains

real(dp) function nen(x, y, z_) result(n)
real(dp), intent(in) :: x, y, z_
real(dp), parameter :: alpha = 6
real(dp) :: r
r = sqrt((x-L/2-L/64)**2+(y-L/2)**2+(z_-L/2)**2)
! This density:
n = -Z*alpha**3/pi**(3._dp/2)*exp(-alpha**2*R**2)
! Corresponds to the potential:
!V = -Z*erf(alpha*R)/R
end function

real(dp) function ne(x, y, z) result(n)
real(dp), intent(in) :: x, y, z
real(dp), parameter :: alpha = 5, Z_ = 1
real(dp) :: r
r = sqrt((x-L/2)**2+(y-L/2)**2+(z-L/2)**2)
n = Z_*alpha**3/pi**(3._dp/2)*exp(-alpha**2*R**2)
end function

end program
