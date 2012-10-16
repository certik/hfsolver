module laplace_assembly
use types, only: dp
implicit none
private
public assemble_2d

real(dp), parameter :: omega = 1.138_dp

contains

real(dp) elemental function f(x, y)
real(dp), intent(in) :: x, y
f = omega**2 * (x**2 + y**2) / 2  ! Harmonic oscillator
end function

subroutine assemble_2d(xin, nodes, elems, ib, xiq, wtq, phihq, dphihq, Am, Bm)
! Assemble on a 2D rectangular uniform mesh
real(dp), intent(in):: xin(:), nodes(:, :), xiq(:), wtq(:, :), &
    phihq(:, :), dphihq(:, :)
integer, intent(in):: elems(:, :), ib(:, :, :)
real(dp), intent(out):: Am(:,:), Bm(:, :)
integer :: Ne, p, e, i, j, iqx, iqy
real(dp), dimension(size(xiq), size(xiq), size(xin), size(xin)) :: &
    phi_v, phi_dx, phi_dy
real(dp) :: x(size(xiq)), y(size(xiq))
real(dp) :: fq(size(xiq), size(xiq)), xp(size(xiq)), yp(size(xiq))
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
                phi_v(iqx, iqy, ax, ay) = phihq(iqx, ax)*phihq(iqy, ay)
                phi_dx(iqx, iqy, ax, ay) = dphihq(iqx, ax)*phihq(iqy, ay)
                phi_dy(iqx, iqy, ax, ay) = phihq(iqx, ax)*dphihq(iqy, ay)
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
    fq = fq * jac_det * wtq
    do by = 1, p+1
    do bx = 1, p+1
        j = ib(bx, by, e)
        if (j==0) cycle
        do ay = 1, p+1
        do ax = 1, p+1
            i = ib(ax, ay, e)
            if (i == 0) cycle
            if (j > i) cycle
            Am(i,j) = Am(i,j) + sum( &
                    (phi_dx(:, :, ax, ay)*phi_dx(:, :, bx, by) &
                    +phi_dy(:, :, ax, ay)*phi_dy(:, :, bx, by)) &
                * jac_det * wtq) / 2
            Am(i,j) = Am(i,j) + sum(fq * &
                phi_v(:, :, ax, ay)*phi_v(:, :, bx, by) &
                * jac_det * wtq)
            Bm(i,j) = Bm(i,j) + sum((phi_v(:, :, ax, ay)*phi_v(:, :, bx, by) &
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

end module




program laplace2d_eig

use types, only: dp
use fe_mesh, only: cartesian_mesh_2d, cartesian_mesh_3d, &
    define_connect_tensor_2d, c2fullc_2d, fe2quad_2d
use laplace_assembly, only: assemble_2d
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
    lam(:), wtq2(:, :), phihq(:, :), dphihq(:, :)
integer, allocatable :: ib(:, :, :), in(:, :, :)
real(dp) :: rmax
integer :: i, j, Nex, Ney

Nex = 1
Ney = 1
p = 3
Nq = 4
rmax = pi/2

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

call assemble_2d(xin, nodes, elems, ib, xiq, wtq2, phihq, dphihq, A, B)
call eigh(A, B, lam, c)
print *, "Eigenvalues:"
do i = 1, min(Nb, 20)
    print *, i, lam(i)
end do
end program
