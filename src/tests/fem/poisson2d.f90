module poisson_assembly
use types, only: dp
use feutils, only: dphih
use linalg, only: inv
use utils, only: assert, stop_error
use constants, only: pi
implicit none
private
public assemble_2d, sol_error

contains

real(dp) elemental function f(x, y)
real(dp), intent(in) :: x, y
f = 2 * pi**2 * sin(pi*x) * sin(pi*y)
end function

real(dp) elemental function exact_sol(x, y)
real(dp), intent(in) :: x, y
exact_sol = sin(pi*x) * sin(pi*y)
end function

subroutine assemble_2d(xin, nodes, elems, ib, xiq, wtq, phihq, Am, rhs)
! Assemble on a 2D rectangular uniform mesh
real(dp), intent(in):: xin(:), nodes(:, :), xiq(:), wtq(:, :), phihq(:, :)
integer, intent(in):: elems(:, :), ib(:, :, :)
real(dp), intent(out):: Am(:,:), rhs(:)
integer :: Ne, Nb, p, e, i, j, iqx, iqy
real(dp) :: dphihq(size(xiq),size(xin))
real(dp), dimension(size(xiq), size(xiq), size(xin), size(xin)) :: &
    phi_v, phi_dx, phi_dy
real(dp) :: x(size(xiq)), y(size(xiq))
real(dp) :: fq(size(xiq), size(xiq)), xp(size(xiq)), yp(size(xiq))
real(dp) :: lx, ly
integer :: ax, ay, bx, by
real(dp) :: jacx, jacy, jac_det

Ne = size(elems, 2)
Nb = maxval(ib)
p = size(xin) - 1
! 1D shape functions
do ax = 1, p+1
    do iqx = 1, size(xiq)
        dphihq(iqx, ax) = dphih(xin, ax, xiq(iqx))
   end do
end do
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
Am=0; rhs=0
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
                * jac_det * wtq)
        end do
        end do
        rhs(j) = rhs(j) + sum(phi_v(:, :, bx, by) * fq)
    end do
    end do
end do
do j = 1, Nb
    do i = 1, j-1
        Am(i, j) = Am(j, i)
    end do
end do
end subroutine

real(dp) pure function det(A)
real(dp), intent(in) :: A(:, :)
det = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)
end function

real(dp) function sol_error(nodes, elems, xiq, wtq, solq) result(r)
real(dp), intent(in) :: nodes(:, :)
integer, intent(in) :: elems(:, :)
real(dp), intent(in) :: xiq(:), wtq(:,:), solq(:, :, :)
real(dp) :: fq(size(xiq), size(xiq)), x(size(xiq)), y(size(xiq))
real(dp) :: jacx, jacy, l(2), jac_det
integer :: e, iqx, iqy, Ne
Ne = size(elems, 2)
r = 0
do e = 1, Ne
    ! l is the diagonal vector:
    l = nodes(:, elems(3, e)) - nodes(:, elems(1, e))
    ! Assume rectangular shape:
    jacx = l(1)/2
    jacy = l(2)/2
    jac_det = abs(jacx*jacy)
    x = (xiq + 1) * jacx + nodes(1, elems(1, e))
    y = (xiq + 1) * jacy + nodes(2, elems(1, e))
    do iqy = 1, size(xiq)
        do iqx = 1, size(xiq)
            fq(iqx, iqy) = exact_sol(x(iqx), y(iqy))
        end do
    end do
    fq = fq-solq(:, :, e)
    r = r + sum(fq*fq * jac_det * wtq)
end do
r = sqrt(r)
end function

real(dp) pure function integrate_3d(f, wtq) result(r)
! Calculates a 3D integral over (-1, 1)^3
real(dp), intent(in) :: f(:, :, :) ! The function f(x,y,z) at quadrature points
real(dp), intent(in) :: wtq(:) ! The 1D quadrature
integer :: i, j, k
r = 0
do i = 1, size(wtq)
    do j = 1, size(wtq)
        do k = 1, size(wtq)
            r = r + f(i, j, k)*wtq(i)*wtq(j)*wtq(k)
        end do
    end do
end do
end function

subroutine append(list, x)
integer, intent(in) :: x
integer, intent(inout) :: list(:)
integer :: i
do i = 1, size(list)
    if (list(i) == 0) then
        list(i) = x
        return
    end if
end do
call stop_error("The list is too short")
end subroutine

integer function get_edge_index(edge_list, v1, v2) result(idx)
integer, intent(inout) :: edge_list(:, :)
integer, intent(in) :: v1, v2
integer :: i, vmin, vmax
vmin = min(v1, v2)
vmax = max(v1, v2)
do i = 1, size(edge_list, 2)
    if (edge_list(1, i) == 0) then
        edge_list(:, i) = [vmin, vmax]
        idx = i
        return
    end if
    if (all(edge_list(:, i) == [vmin, vmax])) then
        idx = i
        return
    end if
end do
idx = -1 ! To fix compiler warning
call stop_error("edge_list is too short")
end function




end module

module poisson2d_code

use types, only: dp
use feutils, only: phih
use fe_mesh, only: cartesian_mesh_2d, cartesian_mesh_3d, &
    define_connect_tensor_2d, c2fullc_2d, fe2quad_2d
use poisson_assembly, only: assemble_2d, sol_error
use feutils, only: get_parent_nodes, get_parent_quad_pts_wts
use linalg, only: solve
use constants, only: pi
implicit none

contains

real(dp) function logg(Nex, Ney, p) result(error)
integer, intent(in) :: p

integer :: Nn, Ne
! nodes(:, i) are the (x,y) coordinates of the i-th mesh node
real(dp), allocatable :: nodes(:, :)
integer, allocatable :: elems(:, :) ! elems(:, i) are nodes of the i-th element
integer :: Nq, Nb
real(dp), allocatable :: xin(:), xiq(:), wtq(:), A(:, :), rhs(:), sol(:), &
        fullsol(:), solq(:, :, :), wtq2(:, :), phihq(:, :)
integer, allocatable :: in(:, :, :), ib(:, :, :)
integer :: i, j
integer, intent(in) :: Nex, Ney

! Test the 3D mesh:
!call cartesian_mesh_3d(3, 10, 2, [0._dp, 0._dp, 0._dp], [1._dp, 2._dp, 3._dp], nodes, elems)

! Test the 2D mesh:
call cartesian_mesh_2d(Nex, Ney, [0.0_dp, 0._dp], [1._dp, 1._dp], nodes, elems)
Nn = size(nodes, 2)
Ne = size(elems, 2)

!print *, "Number of nodes:", Nn
!print *, "Number of elements:", Ne
!print *, "Nodes:"
!do i = 1, Nn
!    print *, i, nodes(:, i)
!end do
!print *, "Elements:"
!do i = 1, Ne
!    print *, i, elems(:, i)
!end do

Nq = p+1
!print *, "Nq =", Nq
!print *, "p =", p
allocate(xin(p+1))
call get_parent_nodes(2, p, xin)
allocate(xiq(Nq), wtq(Nq), wtq2(Nq, Nq))
call get_parent_quad_pts_wts(1, Nq, xiq, wtq)
do j = 1, Nq
    do i = 1, Nq
        wtq2(i, j) = wtq(i)*wtq(j)
    end do
end do
allocate(phihq(size(xiq), size(xin)))
! tabulate parent basis at quadrature points
do j = 1, size(xin)
    do i = 1, size(xiq)
        phihq(i, j) = phih(xin, j, xiq(i))
    end do
end do

call define_connect_tensor_2d(Nex, Ney, p, 1, in)
call define_connect_tensor_2d(Nex, Ney, p, 2, ib)
Nb = maxval(ib)
!print *, "DOFs =", Nb
!print *, in
!print *
!print *, ib
allocate(A(Nb, Nb), rhs(Nb), sol(Nb), fullsol(maxval(in)), solq(Nq, Nq, Ne))

call assemble_2d(xin, nodes, elems, ib, xiq, wtq2, phihq, A, rhs)
sol = solve(A, rhs)
call c2fullc_2d(in, ib, sol, fullsol)
call fe2quad_2d(elems, xin, xiq, phihq, in, fullsol, solq)
error = sol_error(nodes, elems, xiq, wtq2, solq)
!print *, "L2 error:", error

end function

end module




program poisson2d
use types, only: dp
use poisson2d_code, only: logg

real(dp) :: error
integer :: i, N
real(dp) :: t1, t2

N = 4000
call cpu_time(t1)
do i = 1, N
    error = logg(1, 1, 8)
end do
call cpu_time(t2)
print *, "L2 error:", error
print "('Total time: ', f10.8)", (t2-t1) / N
end program
