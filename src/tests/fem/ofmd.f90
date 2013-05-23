module ofmd_utils

use types, only: dp
use feutils, only: phih, dphih
use fe_mesh, only: cartesian_mesh_3d, define_connect_tensor_3d, &
    c2fullc_3d, fe2quad_3d
use poisson3d_assembly, only: assemble_3d, integral, func2quad, func_xyz
use feutils, only: get_parent_nodes, get_parent_quad_pts_wts
use linalg, only: solve
use isolve, only: solve_cg
use utils, only: assert, zeros
use constants, only: pi
implicit none
private
public free_energy

contains

subroutine free_energy(box_dim, Nex, Ney, Nez, p, ibc, frhs, Eh, Nb)
integer, intent(in) :: p, ibc
procedure(func_xyz) :: fexact, frhs
real(dp), intent(in) :: box_dim(3)
real(dp), intent(out) :: Eh
integer, intent(out) :: Nb

integer :: Nn, Ne
! nodes(:, i) are the (x,y) coordinates of the i-th mesh node
real(dp), allocatable :: nodes(:, :)
integer, allocatable :: elems(:, :) ! elems(:, i) are nodes of the i-th element
integer :: Nq
real(dp), allocatable :: xin(:), xiq(:), wtq(:), Ax(:), &
        rhs(:), sol(:), &
        fullsol(:), solq(:, :, :, :), wtq3(:, :, :), phihq(:, :), dphihq(:, :),&
        rhsq(:, :, :, :), exactq(:, :, :, :)
integer, allocatable :: in(:, :, :, :), ib(:, :, :, :), Ap(:), Aj(:)
integer :: i, j, k
integer, intent(in) :: Nex, Ney, Nez
real(dp) :: background

call cartesian_mesh_3d(Nex, Ney, Nez, &
    [0._dp, 0._dp, 0._dp], box_dim, nodes, elems)
Nn = size(nodes, 2)
Ne = size(elems, 2)
Nq = p+1

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
allocate(rhs(Nb), sol(Nb), fullsol(maxval(in)), solq(Nq, Nq, Nq, Ne))
allocate(rhsq(Nq, Nq, Nq, Ne))
allocate(exactq(Nq, Nq, Nq, Ne))

rhsq = func2quad(nodes, elems, xiq, frhs)
! Make the rhsq net neutral (zero integral):
if (ibc == 3) then
    background = integral(nodes, elems, wtq3, rhsq) / &
        (box_dim(1)*box_dim(2)*box_dim(3))
    print *, "Subtracting constant background: ", background
    rhsq = rhsq - background
end if
call assemble_3d(xin, nodes, elems, ib, xiq, wtq3, phihq, dphihq, 4*pi*rhsq, &
    Ap, Aj, Ax, rhs)
print *, "sum(rhs):    ", sum(rhs)
print *, "integral rhs:", integral(nodes, elems, wtq3, rhsq)
print *, "Solving..."
!sol = solve(A, rhs)
sol = solve_cg(Ap, Aj, Ax, rhs, zeros(size(rhs)), 1e-12_dp, 200)
call c2fullc_3d(in, ib, sol, fullsol)
call fe2quad_3d(elems, xin, xiq, phihq, in, fullsol, solq)
if (ibc == 3) solq = solq + (exactq(1, 1, 1, 1) - solq(1, 1, 1, 1))
Eh = integral(nodes, elems, wtq3, solq*rhsq) / 2
print *, "Hartree Energy:", Eh
end subroutine

end module


! ------------------------------------------------------------------------


program ofmd
use types, only: dp
use ofmd_utils, only: free_energy
use constants, only: pi
implicit none
real(dp) :: E
integer :: u, p, DOF, N
real(dp) :: Z, Rcut, Ediff

open(newunit=u, file="H.pseudo", status="old")
read(u, *) Z, N, Rcut, Ediff
p = 4
call free_energy([2._dp, 2._dp, 2._dp], 3, 3, 3, p, 3, rhs, E, DOF)
print *, p, DOF, E
print *, Z, N, Rcut, Ediff

contains

real(dp) function rhs(x, y, z) result(r)
real(dp), intent(in) :: x, y, z
r = 3*pi*exp(sin(pi*x)*sin(pi*y)*sin(pi*z))/4
end function

end program
