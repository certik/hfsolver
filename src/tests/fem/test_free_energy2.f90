module test_free_energy2_utils

use types, only: dp
use feutils, only: phih, dphih
use fe_mesh, only: cartesian_mesh_3d, define_connect_tensor_3d, &
    c2fullc_3d, fe2quad_3d, vtk_save, fe_eval_xyz, line_save, &
    fe2quad_3d_lobatto
use poisson3d_assembly, only: assemble_3d, integral, func2quad, func_xyz, &
    assemble_3d_precalc, assemble_3d_csr, assemble_3d_coo_rhs
use feutils, only: get_parent_nodes, get_parent_quad_pts_wts
!use linalg, only: solve
use isolve, only: solve_cg
use utils, only: assert, zeros, stop_error
use constants, only: pi
use xc, only: xc_pz
use optimize, only: bracket, brent
use ofdft, only: f
use ofdft_fe, only: fe_data, WITH_UMFPACK, free_energy, initialize_fe
use umfpack, only: solve
implicit none
private
public free_energy_min

contains

subroutine free_energy_min(Nelec, L, Nex, Ney, Nez, p, T_au, fnen, fn_pos, &
        Nq, quad_type, &
        Eh, Een, Ts, Exc, Nb)
! Higher level subroutine that does FE precalculation inside it, so it is slow,
! but the interface is simpler.
real(dp), intent(in) :: Nelec, L, T_au
integer, intent(in) :: p, Nq, quad_type, Nex, Ney, Nez
procedure(func_xyz) :: fnen ! (negative) ionic particle density
procedure(func_xyz) :: fn_pos ! (positive) electronic particle density
real(dp), intent(out) :: Eh, Een, Ts, Exc
integer, intent(out) :: Nb

real(dp), allocatable, dimension(:, :, :, :) :: nenq_pos, nq_pos
type(fe_data) :: fed

call initialize_fe(L, Nex, Ney, Nez, p, Nq, quad_type, fed)

allocate(nenq_pos(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(nq_pos(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
nenq_pos = func2quad(fed%nodes, fed%elems, fed%xiq, fnen)
nq_pos = func2quad(fed%nodes, fed%elems, fed%xiq, fn_pos)

call free_energy_min_low_level(Nelec, T_au, nenq_pos, nq_pos, &
        fed, Eh, Een, Ts, Exc)

Nb = fed%Nb
end subroutine

subroutine free_energy_min_low_level(Nelec, T_au, nenq_pos, nq_pos, &
        ! Geometry, finite element specification and internal datastructures
        fed, &
        ! Output arguments
        Eh, Een, Ts, Exc)
! Lower level subroutine that doesn't do any FE precalculation, so it is fast.
real(dp), intent(in) :: Nelec
real(dp), intent(in) :: nenq_pos(:, :, :, :)
real(dp), intent(inout) :: nq_pos(:, :, :, :)
real(dp), intent(in) :: T_au
type(fe_data), intent(in) :: fed
real(dp), intent(out) :: Eh, Een, Ts, Exc

real(dp), allocatable, dimension(:, :, :, :) :: Hpsi, &
    psi, Venq, &
    nenq_neutral
integer :: max_iter
real(dp) :: brent_eps, free_energy_
real(dp) :: psi_norm
real(dp) :: background
real(dp), allocatable :: rhs(:), sol(:), fullsol(:), fullsol2(:)

brent_eps = 1e-3_dp
max_iter = 200

allocate(nenq_neutral(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(Venq(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(Hpsi(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(psi(fed%Nq, fed%Nq, fed%Nq, fed%Ne))


! Calculate Venq
allocate(rhs(fed%Nb), sol(fed%Nb), fullsol(maxval(fed%in)))
allocate(fullsol2(maxval(fed%in)))
background = integral(fed%nodes, fed%elems, fed%wtq3, nenq_pos) / &
    (fed%Lx*fed%Ly*fed%Lz)
print *, "Total (negative) ionic charge: ", background * (fed%Lx*fed%Ly*fed%Lz)
print *, "Subtracting constant background (Q/V): ", background
nenq_neutral = nenq_pos - background
print *, "Assembling RHS..."
call assemble_3d_coo_rhs(fed%Ne, fed%p, 4*pi*nenq_neutral, fed%jac_det, fed%wtq3, &
    fed%ib, fed%phi_v, rhs)
print *, "sum(rhs):    ", sum(rhs)
print *, "integral rhs:", integral(fed%nodes, fed%elems, fed%wtq3, nenq_neutral)
print *, "Solving..."
if (WITH_UMFPACK) then
    call solve(fed%Ap, fed%Aj, fed%Ax, sol, rhs, fed%matd)
else
    sol = solve_cg(fed%Ap, fed%Aj, fed%Ax, rhs, zeros(size(rhs)), 1e-12_dp, 400)
end if
print *, "Converting..."
call c2fullc_3d(fed%in, fed%ib, sol, fullsol)
if (fed%spectral) then
    call fe2quad_3d_lobatto(fed%elems, fed%xiq, fed%in, fullsol, Venq)
else
    call fe2quad_3d(fed%elems, fed%xin, fed%xiq, fed%phihq, fed%in, fullsol, Venq)
end if
print *, "Done"


psi = sqrt(nq_pos)
psi_norm = integral(fed%nodes, fed%elems, fed%wtq3, psi**2)
print *, "Initial norm of psi:", psi_norm
psi = sqrt(Nelec / psi_norm) * psi
psi_norm = integral(fed%nodes, fed%elems, fed%wtq3, psi**2)
print *, "norm of psi:", psi_norm
! This returns H[n] = delta F / delta n, we save it to the Hpsi variable to
! save space:
call free_energy(fed%nodes, fed%elems, fed%in, fed%ib, fed%Nb, fed%Lx, fed%Ly, fed%Lz, fed%xin, fed%xiq, fed%wtq3, T_au, &
    Venq, psi**2, fed%phihq, fed%Ap, fed%Aj, fed%Ax, fed%matd, fed%spectral, &
    fed%phi_v, fed%jac_det, &
    Eh, Een, Ts, Exc, free_energy_, Hpsi=Hpsi)
print *, "Summary of energies [a.u.]:"
print "('    Ts   = ', f14.8)", Ts
print "('    Een  = ', f14.8)", Een
print "('    Eee  = ', f14.8)", Eh
print "('    Exc  = ', f14.8)", Exc
print *, "   ---------------------"
print "('    Etot = ', f14.8, ' a.u.')", free_energy_
call assert(abs(Ts - 10.61904507_dp) < 1e-8_dp)
call assert(abs(Een - (-3.77207363_dp)) < 3e-3_dp)
call assert(abs(Eh - 1.30109486_dp) < 3e-4_dp)
call assert(abs(Exc - (-1.43805889_dp)) < 1e-8_dp)
end subroutine

end module


! ------------------------------------------------------------------------


program test_free_energy2
use types, only: dp
use test_free_energy2_utils, only: free_energy_min
use ofdft, only: read_pseudo
use constants, only: Ha2eV, pi
use utils, only: loadtxt, stop_error
use splines, only: spline3pars, iixmin, poly3, spline3ders
use interp3d, only: trilinear
use feutils, only: quad_gauss
implicit none
real(dp) :: Eh, Een, Ts, Exc
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
!call free_energy_min(L, 4, 4, 4, p, T_au, nen_splines, ne, Eh, Een, Ts, Exc, DOF)
Nq = 20
call free_energy_min(1._dp, L, 4, 4, 4, p, T_au, nen_splines, ne, &
        Nq, quad_gauss, &
        Eh, Een, Ts, Exc, DOF)

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
