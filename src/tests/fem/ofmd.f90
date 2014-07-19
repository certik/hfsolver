program ofmd
use types, only: dp
use ofdft_fe, only: free_energy_min, radial_density_fourier
use ofdft, only: read_pseudo
use constants, only: Ha2eV, pi
use utils, only: loadtxt, stop_error
use splines, only: spline3pars, iixmin, poly3, spline3ders
use interp3d, only: trilinear
use utils, only: linspace
implicit none
real(dp) :: Eh, Een, Ts, Exc, Etot
integer :: p, DOF, u, Nmesh, N, Nx, Ny, Nz, i
real(dp) :: Z, Ediff
real(dp), allocatable :: R(:), V(:), c(:, :)
real(dp), allocatable :: density_en(:), R2(:)
real(dp) :: Rcut, L, T_eV, T_au

p = 6
N = 3
L = 2.997672536043746_dp
T_eV = 0.0862_dp
T_au = T_ev / Ha2eV
Nx = N
Ny = N
Nz = N

call read_pseudo("H.pseudo.gaussian", R, V, Z, Ediff)

Nmesh = 10000
allocate(R2(Nmesh), density_en(Nmesh))
R2 = linspace(0._dp, L/2, Nmesh)
call radial_density_fourier(R, V, L, Z, 32, R2, density_en)

open(newunit=u, file="H.pseudo.density.ft", status="replace")
write(u, "(a)") "# Density. The lines are: r, n(r)"
write(u, *) R2
write(u, *) density_en
close(u)

allocate(c(0:4, size(R2, 1)-1))
call spline3pars(R2, density_en, [2, 2], [0._dp, 0._dp], c)
Rcut = R(size(R))

open(newunit=u, file="H.pseudo.density.ft2", status="replace")
write(u, "(a)") "# Density. The lines are: r, n(r)"
write(u, *) R2
write(u, *) (nen_splines(R2(i), 0._dp, 0._dp), i=1, size(R2))
close(u)

call free_energy_min(L, Nx, Ny, Nz, p, T_au, nen_splines, ne, &
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

contains

real(dp) function nen_splines(x_, y_, z_) result(n)
real(dp), intent(in) :: x_, y_, z_
real(dp) :: r_
integer :: ip
! One atom in the center:
r_ = sqrt((x_+L/64)**2+y_**2+z_**2)
if (r_ >= Rcut) then
    n = 0
    return
end if
!if (r_ <= 1e-4_dp) r_ = 1e-4_dp
ip = 0
ip = iixmin(r_, R2, ip)
n = poly3(r_, c(:, ip))
! FIXME: We need to be using "rho", and flip the sign here:
n = -n
end function

real(dp) function ne(x, y, z) result(n)
real(dp), intent(in) :: x, y, z
real(dp), parameter :: alpha = 0.5_dp, Z_ = 1
real(dp) :: r
r = sqrt(x**2+y**2+z**2)
n = Z_*alpha**3/pi**(3._dp/2)*exp(-alpha**2*R**2)
!n = 1
end function

end program
