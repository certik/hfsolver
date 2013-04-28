module poisson3d_test_util

use types, only: dp
use feutils, only: phih, dphih
use fe_mesh, only: cartesian_mesh_3d, define_connect_tensor_3d, &
    c2fullc_3d, fe2quad_3d
use poisson3d_assembly, only: assemble_3d, integral, func2quad, func_xyz
use feutils, only: get_parent_nodes, get_parent_quad_pts_wts
use linalg, only: solve
use utils, only: assert
implicit none

contains

subroutine test_poisson(box_dim, Nex, Ney, Nez, p, ibc, fexact, frhs, &
        l2_error_eps, hartree_energy_exact, hartree_energy_eps)
integer, intent(in) :: p, ibc
procedure(func_xyz) :: fexact, frhs
real(dp), intent(in) :: l2_error_eps, hartree_energy_exact, hartree_energy_eps
real(dp), intent(in) :: box_dim(3)

integer :: Nn, Ne
! nodes(:, i) are the (x,y) coordinates of the i-th mesh node
real(dp), allocatable :: nodes(:, :)
integer, allocatable :: elems(:, :) ! elems(:, i) are nodes of the i-th element
integer :: Nq, Nb
real(dp), allocatable :: xin(:), xiq(:), wtq(:), A(:, :), rhs(:), sol(:), &
        fullsol(:), solq(:, :, :, :), wtq3(:, :, :), phihq(:, :), dphihq(:, :),&
        rhsq(:, :, :, :), exactq(:, :, :, :)
integer, allocatable :: in(:, :, :, :), ib(:, :, :, :)
integer :: i, j, k
integer, intent(in) :: Nex, Ney, Nez
real(dp) :: hartree_energy, l2_error, hartree_energy_error

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
allocate(A(Nb, Nb), rhs(Nb), sol(Nb), fullsol(maxval(in)), solq(Nq, Nq, Nq, Ne))
allocate(rhsq(Nq, Nq, Nq, Ne))
allocate(exactq(Nq, Nq, Nq, Ne))

exactq = func2quad(nodes, elems, xiq, fexact)
rhsq = func2quad(nodes, elems, xiq, frhs)
call assemble_3d(xin, nodes, elems, ib, xiq, wtq3, phihq, dphihq, rhsq, A, rhs)
print *, "Solving..."
sol = solve(A, rhs)
call c2fullc_3d(in, ib, sol, fullsol)
call fe2quad_3d(elems, xin, xiq, phihq, in, fullsol, solq)
l2_error = sqrt(integral(nodes, elems, wtq3, (solq-exactq)**2))
hartree_energy = integral(nodes, elems, wtq3, solq*rhsq)
hartree_energy_error = abs(hartree_energy-hartree_energy_exact)
print "(' L2 Error:',es10.2)", l2_error
print *, "Hartree Energy (calculated):", hartree_energy
print *, "Hartree Energy (exact):     ", hartree_energy_exact
print "(' Hartree Energy (error):',es10.2)", hartree_energy_error
call assert(l2_error < l2_error_eps)
call assert(hartree_energy_error < hartree_energy_eps)
end subroutine

end module


! ------------------------------------------------------------------------


program test_poisson3d
use types, only: dp
use poisson3d_test_util, only: test_poisson
use constants, only: pi
implicit none

call test_poisson([1._dp, 1._dp, 1._dp], 1, 1, 1, 8, 2, sol, rhs, &
    1e-7_dp, 3*pi**2/8, 1e-11_dp)
call test_poisson([1._dp, 1._dp, 1._dp], 2, 3, 5, 5, 2, sol, rhs, &
    1e-5_dp, 3*pi**2/8, 1e-7_dp)

call test_poisson([2._dp, 3._dp, 5._dp], 2, 2, 2, 8, 2, sol, rhs, &
    1e-2_dp, 45*pi**2/4, 1e-3_dp)

call test_poisson([2._dp, 3._dp, 5._dp], 2, 3, 5, 6, 2, sol, rhs, &
    1e-4_dp, 45*pi**2/4, 1e-5_dp)

call test_poisson([2._dp, 2._dp, 2._dp], 2, 3, 5, 6, 3, sol, rhs, &
    8e-2_dp, 3*pi**2, 5e-7_dp)

call test_poisson([2._dp, 2._dp, 2._dp], 2, 3, 5, 6, 3, sol2, rhs2, &
    1e-1_dp, 3*pi**2, 5e-6_dp)

contains

real(dp) function rhs(x, y, z) result(r)
real(dp), intent(in) :: x, y, z
r = 3 * pi**2 * sin(pi*x) * sin(pi*y) * sin(pi*z)
end function

real(dp) function sol(x, y, z) result(r)
real(dp), intent(in) :: x, y, z
r = sin(pi*x) * sin(pi*y) * sin(pi*z)
end function

real(dp) function rhs2(x, y, z) result(r)
real(dp), intent(in) :: x, y, z
r = 3 * pi**2 * sin(pi*(x+0.5_dp)) * sin(pi*(y+0.5_dp)) * sin(pi*(z+0.5_dp))
end function

real(dp) function sol2(x, y, z) result(r)
real(dp), intent(in) :: x, y, z
r = sin(pi*(x+0.5_dp)) * sin(pi*(y+0.5_dp)) * sin(pi*(z+0.5_dp))
end function

end program
