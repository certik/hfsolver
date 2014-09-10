program test_free_energy4

! Test for 4 Gaussian charges.

! This test uses FE and produces the same result as test_free_energy_fft2

use types, only: dp
use ofdft_fe, only: free_energy2
use constants, only: Ha2eV, pi
use utils, only: loadtxt, assert
use splines, only: spline3pars, iixmin, poly3
use interp3d, only: trilinear
use feutils, only: quad_lobatto
use md, only: positions_fcc
implicit none
real(dp) :: Eee, Een, Ts, Exc, Etot
integer :: p, DOF, Nq
real(dp) :: Rcut, L, T_eV, T_au
integer, parameter :: natom = 4
real(dp) :: X(3, natom)

Rcut = 0.3_dp
p = 8
L = 2
T_eV = 0.0862_dp
T_au = T_ev / Ha2eV
Nq = 9
call positions_fcc(X, L)
call free_energy2(real(natom, dp), L, 8, 8, 8, p, T_au, nen, ne, &
        Nq, quad_lobatto, &
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

real(dp) function ne(x, y, z) result(n)
real(dp), intent(in) :: x, y, z
real(dp), parameter :: alpha = 5, Z_ = 4
real(dp) :: r
r = sqrt((x-L/2)**2+(y-L/2)**2+(z-L/2)**2)
n = Z_*alpha**3/pi**(3._dp/2)*exp(-alpha**2*R**2)
end function

end program
