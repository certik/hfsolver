module laplace_assembly
use types, only: dp
use feutils, only: phih, dphih
use linalg, only: inv
use utils, only: assert
implicit none
private
public assemble_2d, define_connect

contains

subroutine assemble_2d(xin, nodes, elems, ib, xiq, wtq, Am, Bm)
real(dp), intent(in):: xin(:), nodes(:, :), xiq(:), wtq(:)
integer, intent(in):: elems(:, :), ib(:, :, :)
real(dp), intent(out):: Am(:,:), Bm(:, :)
integer :: Ne, Nb, p, e, i, j, iqx, iqy
real(dp) :: phihq(size(xiq),size(xin)), dphihq(size(xiq),size(xin))
real(dp), dimension(size(xiq), size(xiq), size(xin), size(xin)) :: &
    phi_v, phi_dx, phi_dy
real(dp) :: hq(size(xiq), size(xiq))
real(dp) :: l(2)
integer :: ax, ay, bx, by
real(dp) :: jacx, jacy, jac_det

Ne = size(elems, 2)
Nb = maxval(ib)
p = size(xin) - 1
! 1D shape functions
do ax = 1, p+1
    do iqx = 1, size(xiq)
        phihq(iqx, ax) = phih(xin,ax,xiq(iqx))
        dphihq(iqx, ax) = dphih(xin,ax,xiq(iqx))
   end do
end do
! 2D shape functions
do ax = 1, p+1
    do ay = 1, p+1
        do iqx = 1, size(xiq)
            do iqy = 1, size(xiq)
                phi_v(iqx, iqy, ax, ay) = phihq(iqx, ax)*phihq(iqy,ay)
                phi_dx(iqx, iqy, ax, ay) = dphihq(iqx, ax)*phihq(iqy,ay)
                phi_dy(iqx, iqy, ax, ay) = phihq(iqx, ax)*dphihq(iqy,ay)
            end do
        end do
   end do
end do
Am=0; Bm=0
do e = 1, Ne
    ! l is the diagonal vector:
    l = nodes(:, elems(1, e)) - nodes(:, elems(3, e))
    ! Assume rectangular shape:
    jacx = l(1)/2
    jacy = l(2)/2
    jac_det = abs(jacx*jacy)
    do bx = 1, p+1
    do by = 1, p+1
        j = ib(bx, by, e)
        if (j==0) cycle
        do ax = 1, p+1
        do ay = 1, p+1
            i = ib(ax, ay, e)
            if (i == 0) cycle
            if (j > i) cycle
            hq = phi_dx(:, :, ax, ay)*phi_dx(:, :, bx, by) / (jacx*jacy) &
                +phi_dy(:, :, ax, ay)*phi_dy(:, :, bx, by) / (jacx*jacy)
            hq = hq / 2
            Am(i,j) = Am(i,j) + integrate_2d(hq * jac_det, wtq)
            hq = phi_v(:, :, ax, ay) * phi_v(:, :, bx, by)
            Bm(i,j) = Bm(i,j) + integrate_2d(hq * jac_det, wtq)
        end do
        end do
    end do
    end do
end do
do j = 1, Nb
    do i = 1, j-1
        Am(i, j) = Am(j, i)
        Bm(i, j) = Bm(j, i)
    end do
end do
end subroutine

real(dp) pure function det(A)
real(dp), intent(in) :: A(:, :)
det = A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)
end function

real(dp) pure function integrate_2d(f, wtq) result(r)
! Calculates a 2D integral over (-1, 1)^2
real(dp), intent(in) :: f(:, :) ! The function f(x, y) at quadrature points
real(dp), intent(in) :: wtq(:) ! The 1D quadrature
integer :: i, j
r = 0
do i = 1, size(wtq)
    do j = 1, size(wtq)
        r = r + f(i, j)*wtq(i)*wtq(j)
    end do
end do
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

subroutine define_connect(elems, p, ib)
integer, intent(in):: elems(:, :), p
integer, allocatable, intent(out) :: ib(:, :, :)
integer :: i, j, Ne, idx
Ne = size(elems, 2)
allocate(ib(p+1, p+1, Ne))
call assert(Ne == 1)
ib = 0
idx = 1
do i = 2, p
    do j = 2, p
        ib(i, j, 1) = idx
        idx = idx + 1
    end do
end do
end subroutine

end module




program laplace2d_eig

use types, only: dp
use fe_mesh, only: cartesian_mesh_2d, cartesian_mesh_3d
use laplace_assembly, only: assemble_2d, define_connect
use feutils, only: get_parent_nodes, get_parent_quad_pts_wts
use linalg, only: eigh
use constants, only: pi
implicit none

integer :: Nn, Ne
! nodes(:, i) are the (x,y) coordinates of the i-th mesh node
real(dp), allocatable :: nodes(:, :)
integer, allocatable :: elems(:, :) ! elems(:, i) are nodes of the i-th element
integer :: Nq, p, Nb
real(dp), allocatable :: xin(:), xiq(:), wtq(:), A(:, :), B(:, :), c(:, :), &
    lam(:)
integer, allocatable :: ib(:, :, :)
integer :: i


! Test the 3D mesh:
call cartesian_mesh_3d(3, 10, 2, [0._dp, 0._dp, 0._dp], [1._dp, 2._dp, 3._dp], nodes, elems)

! Test the 2D mesh:
call cartesian_mesh_2d(1, 1, [0._dp, 0._dp], [pi, pi], nodes, elems)
Nn = size(nodes, 2)
Ne = size(elems, 2)

print *, "Number of nodes:", Nn
print *, "Number of elements:", Ne
!print *, "Nodes:"
!do i = 1, Nn
!    print *, i, nodes(:, i)
!end do
!print *, "Elements:"
!do i = 1, Ne
!    print *, i, elems(:, i)
!end do

Nq = 11
p = 10
allocate(xin(p+1))
call get_parent_nodes(2, p, xin)
allocate(xiq(Nq), wtq(Nq))
call get_parent_quad_pts_wts(1, Nq, xiq, wtq)

call define_connect(elems, p, ib)
Nb = maxval(ib)
print *, "DOFs =", Nb
allocate(A(Nb, Nb), B(Nb, Nb), c(Nb, Nb), lam(Nb))

call assemble_2d(xin, nodes, elems, ib, xiq, wtq, A, B)
!print *, "A:"
!print "(4f10.6)", A
!print *, "B:"
!print "(4f10.6)", B
call eigh(A, B, lam, c)
print *, "Eigenvalues:"
do i = 1, Nb
    print *, i, lam(i)
end do
end program
