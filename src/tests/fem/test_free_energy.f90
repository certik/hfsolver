module test_free_energy_utils

use types, only: dp
use feutils, only: phih, dphih
use fe_mesh, only: cartesian_mesh_3d, define_connect_tensor_3d, &
    c2fullc_3d, fe2quad_3d, vtk_save, fe_eval_xyz, line_save
use poisson3d_assembly, only: assemble_3d, integral, func2quad, func_xyz
use feutils, only: get_parent_nodes, get_parent_quad_pts_wts
use linalg, only: solve
use isolve, only: solve_cg
use utils, only: assert, zeros, stop_error
use constants, only: pi
use xc, only: xc_pz
use ofdft, only: f
implicit none
private
public free_energy

contains

subroutine free_energy(L, Nex, Ney, Nez, p, T_au, fnen, fn_pos, Eh, Een, Ts, &
        Exc, Nb)
integer, intent(in) :: p
procedure(func_xyz) :: fnen ! (negative) ionic particle density
procedure(func_xyz) :: fn_pos ! (positive) electronic particle density
real(dp), intent(in) :: L, T_au
real(dp), intent(out) :: Eh, Een, Ts, Exc
integer, intent(out) :: Nb

integer :: Nn, Ne, ibc
! nodes(:, i) are the (x,y) coordinates of the i-th mesh node
real(dp), allocatable :: nodes(:, :)
integer, allocatable :: elems(:, :) ! elems(:, i) are nodes of the i-th element
integer :: Nq
real(dp), allocatable :: xin(:), xiq(:), wtq(:), Ax(:), &
        rhs(:), sol(:), &
        fullsol(:), Vhq(:, :, :, :), wtq3(:, :, :), phihq(:, :), dphihq(:, :),&
        nenq(:, :, :, :), &
        Venq(:, :, :, :), y(:, :, :, :), F0(:, :, :, :), &
        exc_density(:, :, :, :), nq_pos(:, :, :, :), nq_neutral(:, :, :, :)
integer, allocatable :: in(:, :, :, :), ib(:, :, :, :), Ap(:), Aj(:)
integer :: i, j, k, m
integer, intent(in) :: Nex, Ney, Nez
real(dp) :: background
real(dp) :: Lx, Ly, Lz
real(dp) :: beta, tmp

ibc = 3 ! Periodic boundary condition

Lx = L
Ly = L
Lz = L

call cartesian_mesh_3d(Nex, Ney, Nez, &
    [-Lx/2, -Ly/2, -Lz/2], [Lx/2, Ly/2, Lz/2], nodes, elems)
Nn = size(nodes, 2)
Ne = size(elems, 2)
Nq = 20

print *, "Number of nodes:", Nn
print *, "Number of elements:", Ne
print *, "Nq =", Nq
print *, "p =", p
allocate(xin(p+1))
call get_parent_nodes(2, p, xin)
allocate(xiq(Nq), wtq(Nq), wtq3(Nq, Nq, Nq))
call get_parent_quad_pts_wts(1, Nq, xiq, wtq)
forall(i=1:Nq, j=1:Nq, k=1:Nq) wtq3(i, j, k) = wtq(i)*wtq(j)*wtq(k)
allocate(phihq(size(xiq), size(xin)))
allocate(dphihq(size(xiq), size(xin)))
! Tabulate parent basis at quadrature points
forall(i=1:size(xiq), j=1:size(xin))  phihq(i, j) =  phih(xin, j, xiq(i))
forall(i=1:size(xiq), j=1:size(xin)) dphihq(i, j) = dphih(xin, j, xiq(i))

call define_connect_tensor_3d(Nex, Ney, Nez, p, 1, in)
call define_connect_tensor_3d(Nex, Ney, Nez, p, ibc, ib)
Nb = maxval(ib)
print *, "DOFs =", Nb
allocate(rhs(Nb), sol(Nb), fullsol(maxval(in)), Vhq(Nq, Nq, Nq, Ne))
allocate(Venq(Nq, Nq, Nq, Ne))
allocate(nenq(Nq, Nq, Nq, Ne))
allocate(y(Nq, Nq, Nq, Ne))
allocate(F0(Nq, Nq, Nq, Ne))
allocate(exc_density(Nq, Nq, Nq, Ne))
allocate(nq_pos(Nq, Nq, Nq, Ne))
allocate(nq_neutral(Nq, Nq, Nq, Ne))

nenq = func2quad(nodes, elems, xiq, fnen)
nq_pos = func2quad(nodes, elems, xiq, fn_pos)
! Make the charge density net neutral (zero integral):
background = integral(nodes, elems, wtq3, nq_pos) / (Lx*Ly*Lz)
print *, "Total (positive) electronic charge: ", background * (Lx*Ly*Lz)
print *, "Subtracting constant background (Q/V): ", background
nq_neutral = nq_pos - background
call assemble_3d(xin, nodes, elems, ib, xiq, wtq3, phihq, dphihq, &
    4*pi*nq_neutral, Ap, Aj, Ax, rhs)
print *, "sum(rhs):    ", sum(rhs)
print *, "integral rhs:", integral(nodes, elems, wtq3, nq_neutral)
print *, "Solving..."
!sol = solve(A, rhs)
sol = solve_cg(Ap, Aj, Ax, rhs, zeros(size(rhs)), 1e-12_dp, 400)
call c2fullc_3d(in, ib, sol, fullsol)
call fe2quad_3d(elems, xin, xiq, phihq, in, fullsol, Vhq)

background = integral(nodes, elems, wtq3, nenq) / (Lx*Ly*Lz)
print *, "Total (negative) ionic charge: ", background * (Lx*Ly*Lz)
print *, "Subtracting constant background (Q/V): ", background
nenq = nenq - background
call assemble_3d(xin, nodes, elems, ib, xiq, wtq3, phihq, dphihq, &
    4*pi*nenq, Ap, Aj, Ax, rhs)
print *, "sum(rhs):    ", sum(rhs)
print *, "integral rhs:", integral(nodes, elems, wtq3, nenq)
print *, "Solving..."
sol = solve_cg(Ap, Aj, Ax, rhs, zeros(size(rhs)), 1e-12_dp, 400)
call c2fullc_3d(in, ib, sol, fullsol)
call fe2quad_3d(elems, xin, xiq, phihq, in, fullsol, Venq)
print *, "Saving Ven to VTK"
call vtk_save("Venq.vtk", Nex, Ney, Nez, nodes, elems, xiq, Venq)
print *, "Saving values of Ven on a line"
call line_save("Venq_line.txt", xin, nodes, elems, in, fullsol, &
    [0._dp, 0._dp, 0._dp], [1._dp, 0._dp, 0._dp], 500)
print *, "    Done."


! Hartree energy
Eh = integral(nodes, elems, wtq3, Vhq*nq_neutral) / 2
! Electron-nucleus energy
Een = integral(nodes, elems, wtq3, Venq*nq_neutral)
! Kinetic energy using Perrot parametrization
beta = 1/T_au
! The density must be positive, the f(y) fails for negative "y". Thus we use
! nq_pos.
y = pi**2 / sqrt(2._dp) * beta**(3._dp/2) * nq_pos
if (any(y < 0)) call stop_error("Density must be positive")
F0 = nq_pos / beta * f(y)
Ts = integral(nodes, elems, wtq3, F0)
! Exchange and correlation energy
do m = 1, Ne
do k = 1, Nq
do j = 1, Nq
do i = 1, Nq
    call xc_pz(nq_pos(i, j, k, m), exc_density(i, j, k, m), tmp)
end do
end do
end do
end do
Exc = integral(nodes, elems, wtq3, exc_density * nq_pos)
end subroutine

end module


! ------------------------------------------------------------------------


program test_free_energy
use types, only: dp
use test_free_energy_utils, only: free_energy
use constants, only: Ha2eV, pi
use utils, only: loadtxt, assert
use splines, only: spline3pars, iixmin, poly3
use interp3d, only: trilinear
implicit none
real(dp) :: Eh, Een, Ts, Exc, Etot
integer :: p, DOF
real(dp) :: Z
real(dp) :: Rcut, L, T_eV, T_au

Z = 1
Rcut = 0.3_dp
p = 4
L = 2
T_eV = 0.0862_dp
T_au = T_ev / Ha2eV
call free_energy(L, 8, 8, 8, p, T_au, nen, ne, Eh, Een, Ts, Exc, DOF)
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

! The reference answers are converged to at least 1e-5:
print *, abs(Ts  - (+10.61905))
call assert(abs(Ts  - (+10.61905)) < 1e-4)
print *, abs(Een - (- 3.80769))
call assert(abs(Een - (- 3.80769)) < 1e-4)
print *, abs(Eh  - (+ 1.30109))
call assert(abs(Eh  - (+ 1.30109)) < 1e-4)
print *, abs(Exc  - (- 1.43806))
call assert(abs(Exc  - (- 1.43806)) < 1e-4)

contains

real(dp) function nen(x, y, z_) result(n)
real(dp), intent(in) :: x, y, z_
real(dp), parameter :: alpha = 12
real(dp) :: r
r = sqrt(x**2+y**2+z_**2)
! This density:
n = -Z*alpha**3/pi**(3._dp/2)*exp(-alpha**2*R**2)
! Corresponds to the potential:
!V = -Z*erf(alpha*R)/R
end function

real(dp) function ne(x, y, z) result(n)
real(dp), intent(in) :: x, y, z
real(dp), parameter :: alpha = 5, Z_ = 1
real(dp) :: r
r = sqrt(x**2+y**2+z**2)
n = Z_*alpha**3/pi**(3._dp/2)*exp(-alpha**2*R**2)
end function

end program
