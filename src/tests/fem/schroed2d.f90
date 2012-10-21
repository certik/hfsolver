module schroed2d_assembly
use types, only: dp
use sorting, only: sort
implicit none
private
public assemble_2d, Z, exact_energies

real(dp), parameter :: Z = 2._dp

contains

real(dp) elemental function f(x, y)
real(dp), intent(in) :: x, y
f = -Z/sqrt(x**2 + y**2)
end function

subroutine assemble_2d(xin, nodes, elems, ib, xiq, wtq, phihq, dphihq, Am, Bm)
! Assemble on a 2D rectangular uniform mesh
real(dp), intent(in):: xin(:), nodes(:, :), xiq(:), wtq(:, :), &
    phihq(:, :), dphihq(:, :)
integer, intent(in):: elems(:, :), ib(:, :, :)
real(dp), intent(out):: Am(:,:), Bm(:, :)
real(dp), dimension(size(xiq), size(xiq), &
    size(xin), size(xin)) :: phi_v, phi_dx, phi_dy
real(dp), dimension(size(xiq), size(xiq)) :: fq
real(dp), dimension(size(xiq)) :: x, y, xp, yp
integer :: Ne, p, e, i, j, iqx, iqy
real(dp) :: lx, ly
integer :: ax, ay, bx, by
real(dp) :: jacx, jacy, jac_det

Ne = size(elems, 2)
p = size(xin) - 1
! 2D shape functions
do ay = 1, p+1
do ax = 1, p+1
    do iqy = 1, size(xiq)
    do iqx = 1, size(xiq)
        phi_v (iqx, iqy, ax, ay) =  phihq(iqx, ax) *  phihq(iqy, ay)
        phi_dx(iqx, iqy, ax, ay) = dphihq(iqx, ax) *  phihq(iqy, ay)
        phi_dy(iqx, iqy, ax, ay) =  phihq(iqx, ax) * dphihq(iqy, ay)
    end do
    end do
end do
end do
Am=0; Bm=0
! Precalculate as much as possible:
lx = nodes(1, elems(3, 1)) - nodes(1, elems(1, 1)) ! Element sizes
ly = nodes(2, elems(3, 1)) - nodes(2, elems(1, 1))
jacx = lx/2
jacy = ly/2
jac_det = abs(jacx*jacy)
xp = (xiq + 1) * jacx
yp = (xiq + 1) * jacy
phi_dx = phi_dx / jacx
phi_dy = phi_dy / jacy
do e = 1, Ne
    x = xp + nodes(1, elems(1, e))
    y = yp + nodes(2, elems(1, e))
    do iqy = 1, size(xiq)
    do iqx = 1, size(xiq)
        fq(iqx, iqy) = f(x(iqx), y(iqy))
    end do
    end do
    do by = 1, p+1
    do bx = 1, p+1
        j = ib(bx, by, e)
        if (j==0) cycle
        do ay = 1, p+1
        do ax = 1, p+1
            i = ib(ax, ay, e)
            if (i == 0) cycle
            if (j > i) cycle
            Am(i,j) = Am(i,j) + sum(( &
                phi_dx(:, :, ax, ay)*phi_dx(:, :, bx, by) + &
                phi_dy(:, :, ax, ay)*phi_dy(:, :, bx, by)) &
                * jac_det * wtq) / 2
            Am(i,j) = Am(i,j) + sum(fq * &
                phi_v(:, :, ax, ay)*phi_v(:, :, bx, by) &
                * jac_det * wtq)
            Bm(i,j) = Bm(i,j) + sum(( &
                phi_v(:, :, ax, ay) * phi_v(:, :, bx, by) &
                * jac_det * wtq))
        end do
        end do
    end do
    end do
end do
do j = 1, size(Am, 2)
    do i = 1, j-1
        Am(i, j) = Am(j, i)
        Bm(i, j) = Bm(j, i)
    end do
end do
end subroutine

integer pure function calc_size(nmax) result(s)
integer, intent(in) :: nmax
integer :: n, m
s = 0
do n = 1, nmax
    do m = -n+1, n-1
        s = s + 1
    end do
end do
end function

subroutine exact_energies(Z, E)
real(dp), intent(in) :: Z
real(dp), intent(out) :: E(:)
real(dp) :: tmp(calc_size(size(E)))
integer :: n, m, idx
idx = 0
do n = 1, size(E)
    do m = -n+1, n-1
        idx = idx + 1
        tmp(idx) = -Z**2 / (2*(n-1._dp/2)**2)
    end do
end do
E = tmp(:size(E))
end subroutine

end module


!------------------------------------------------------------------------------

program schroed2d

use types, only: dp
use fe_mesh, only: cartesian_mesh_2d, define_connect_tensor_2d
use schroed2d_assembly, only: assemble_2d, Z, exact_energies
use feutils, only: get_parent_nodes, get_parent_quad_pts_wts, phih, dphih
use linalg, only: eigh
use constants, only: pi
implicit none

integer :: Nn, Ne
! nodes(:, i) are the (x,y) coordinates of the i-th mesh node
real(dp), allocatable :: nodes(:, :)
integer, allocatable :: elems(:, :) ! elems(:, i) are nodes of the i-th element
integer :: Nq, p, Nb
real(dp), allocatable :: xin(:), xiq(:), wtq(:), A(:, :), B(:, :), c(:, :), &
    lam(:), wtq2(:, :), phihq(:, :), dphihq(:, :), E_exact(:)
integer, allocatable :: ib(:, :, :), in(:, :, :)
real(dp) :: rmax
integer :: i, j, Nex, Ney

Nex = 2
Ney = 2
p = 20
Nq = p+1
rmax = 15  ! The size of the box in atomic units

call cartesian_mesh_2d(Nex, Ney, [-rmax, -rmax], [rmax, rmax], nodes, elems)
Nn = size(nodes, 2)
Ne = size(elems, 2)

print *, "Number of nodes:", Nn
print *, "Number of elements:", Ne

allocate(xin(p+1))
call get_parent_nodes(2, p, xin)
allocate(xiq(Nq), wtq(Nq), wtq2(Nq, Nq))
call get_parent_quad_pts_wts(1, Nq, xiq, wtq)
forall(i=1:Nq, j=1:Nq) wtq2(i, j) = wtq(i)*wtq(j)
allocate(phihq(size(xiq), size(xin)))
allocate(dphihq(size(xiq), size(xin)))
! Tabulate parent basis at quadrature points
forall(i=1:size(xiq), j=1:size(xin))  phihq(i, j) =  phih(xin, j, xiq(i))
forall(i=1:size(xiq), j=1:size(xin)) dphihq(i, j) = dphih(xin, j, xiq(i))

call define_connect_tensor_2d(Nex, Ney, p, 1, in)
call define_connect_tensor_2d(Nex, Ney, p, 2, ib)
Nb = maxval(ib)
print *, "p =", p
print *, "DOFs =", Nb
allocate(A(Nb, Nb), B(Nb, Nb), c(Nb, Nb), lam(Nb))

print *, "Assembling..."
call assemble_2d(xin, nodes, elems, ib, xiq, wtq2, phihq, dphihq, A, B)
print *, "Solving..."
call eigh(A, B, lam, c)
print *, "Eigenvalues:"
allocate(E_exact(20))
call exact_energies(Z, E_exact)
do i = 1, min(Nb, 20)
    print "(i4, f12.6, f12.6, es12.2)", i, lam(i), E_exact(i), rel(lam(i), E_exact(i))
end do

contains

real(dp) pure function rel(a, b)
real(dp), intent(in) :: a, b
rel = abs(a-b) / max(abs(a), abs(b))
end function

end program
