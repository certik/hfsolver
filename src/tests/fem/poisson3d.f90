module poisson3d_assembly
use types, only: dp
use linalg, only: inv
use utils, only: assert, stop_error
use constants, only: pi
implicit none
private
public assemble_3d, sol_error, integral

contains

real(dp) elemental function f(x, y, z)
real(dp), intent(in) :: x, y, z
f = 3 * pi**2 * sin(pi*x) * sin(pi*y) * sin(pi*z)
end function

real(dp) elemental function exact_sol(x, y, z)
real(dp), intent(in) :: x, y, z
exact_sol = sin(pi*x) * sin(pi*y) * sin(pi*z)
end function

subroutine assemble_3d(xin, nodes, elems, ib, xiq, wtq, phihq, dphihq, Am, rhs)
! Assemble on a 2D rectangular uniform mesh
real(dp), intent(in):: xin(:), nodes(:, :), xiq(:), wtq(:, :, :), &
    phihq(:, :), dphihq(:, :)
integer, intent(in):: elems(:, :), ib(:, :, :, :)
real(dp), intent(out):: Am(:,:), rhs(:)
integer :: Ne, p, e, i, j, iqx, iqy, iqz
real(dp), dimension(size(xiq), size(xiq), size(xiq), &
    size(xin), size(xin), size(xin)) :: phi_v, phi_dx, phi_dy, phi_dz
real(dp), dimension(size(xiq), size(xiq), size(xiq)) :: fq
real(dp), dimension(size(xiq)) :: x, y, z, xp, yp, zp
real(dp) :: lx, ly, lz
integer :: ax, ay, az, bx, by, bz
real(dp) :: jacx, jacy, jacz, jac_det

Ne = size(elems, 2)
p = size(xin) - 1
! 3D shape functions
do az = 1, p+1
do ay = 1, p+1
do ax = 1, p+1
    do iqz = 1, size(xiq)
    do iqy = 1, size(xiq)
    do iqx = 1, size(xiq)
        phi_v (iqx, iqy, iqz, ax, ay, az) = &
             phihq(iqx, ax) *  phihq(iqy, ay) *  phihq(iqz, az)
        phi_dx(iqx, iqy, iqz, ax, ay, az) = &
            dphihq(iqx, ax) *  phihq(iqy, ay) *  phihq(iqz, az)
        phi_dy(iqx, iqy, iqz, ax, ay, az) = &
             phihq(iqx, ax) * dphihq(iqy, ay) *  phihq(iqz, az)
        phi_dz(iqx, iqy, iqz, ax, ay, az) = &
             phihq(iqx, ax) *  phihq(iqy, ay) * dphihq(iqz, az)
    end do
    end do
    end do
end do
end do
end do
Am=0; rhs=0
! Precalculate as much as possible:
lx = nodes(1, elems(7, 1)) - nodes(1, elems(1, 1)) ! Element sizes
ly = nodes(2, elems(7, 1)) - nodes(2, elems(1, 1))
lz = nodes(3, elems(7, 1)) - nodes(3, elems(1, 1))
jacx = lx/2
jacy = ly/2
jacz = lz/2
jac_det = abs(jacx*jacy*jacz)
xp = (xiq + 1) * jacx
yp = (xiq + 1) * jacy
zp = (xiq + 1) * jacz
phi_dx = phi_dx / jacx
phi_dy = phi_dy / jacy
phi_dz = phi_dz / jacz
do e = 1, Ne
    x = xp + nodes(1, elems(1, e))
    y = yp + nodes(2, elems(1, e))
    z = zp + nodes(3, elems(1, e))
    do iqz = 1, size(xiq)
    do iqy = 1, size(xiq)
    do iqx = 1, size(xiq)
        fq(iqx, iqy, iqz) = f(x(iqx), y(iqy), z(iqz))
    end do
    end do
    end do
    fq = fq * jac_det * wtq
    do bz = 1, p+1
    do by = 1, p+1
    do bx = 1, p+1
        j = ib(bx, by, bz, e)
        if (j==0) cycle
        do az = 1, p+1
        do ay = 1, p+1
        do ax = 1, p+1
            i = ib(ax, ay, az, e)
            if (i == 0) cycle
            if (j > i) cycle
            Am(i,j) = Am(i,j) + sum( &
                (phi_dx(:, :, :, ax, ay, az)*phi_dx(:, :, :, bx, by, bz) &
                +phi_dy(:, :, :, ax, ay, az)*phi_dy(:, :, :, bx, by, bz) &
                +phi_dz(:, :, :, ax, ay, az)*phi_dz(:, :, :, bx, by, bz)) &
                    * jac_det * wtq)
        end do
        end do
        end do
        rhs(j) = rhs(j) + sum(phi_v(:, :, :, bx, by, bz) * fq)
    end do
    end do
    end do
end do
do j = 1, size(Am, 2)
    do i = 1, j-1
        Am(i, j) = Am(j, i)
    end do
end do
end subroutine

real(dp) function sol_error(nodes, elems, xiq, wtq, solq) result(r)
real(dp), intent(in) :: nodes(:, :)
integer, intent(in) :: elems(:, :)
real(dp), intent(in) :: xiq(:), wtq(:, :, :), solq(:, :, :, :)
real(dp) :: fq(size(xiq), size(xiq), size(xiq)), x(size(xiq)), y(size(xiq)), &
    z(size(xiq))
real(dp) :: jacx, jacy, jacz, l(3), jac_det
integer :: e, iqx, iqy, iqz, Ne
Ne = size(elems, 2)
r = 0
do e = 1, Ne
    ! l is the diagonal vector:
    l = nodes(:, elems(7, e)) - nodes(:, elems(1, e))
    ! Assume rectangular shape:
    jacx = l(1)/2
    jacy = l(2)/2
    jacz = l(3)/2
    jac_det = abs(jacx*jacy*jacz)
    x = (xiq + 1) * jacx + nodes(1, elems(1, e))
    y = (xiq + 1) * jacy + nodes(2, elems(1, e))
    z = (xiq + 1) * jacz + nodes(3, elems(1, e))
    do iqz = 1, size(xiq)
    do iqy = 1, size(xiq)
    do iqx = 1, size(xiq)
            fq(iqx, iqy, iqz) = exact_sol(x(iqx), y(iqy), z(iqz))
    end do
    end do
    end do
    fq = fq-solq(:, :, :, e)
    r = r + sum(fq*fq * jac_det * wtq)
end do
r = sqrt(r)
end function

real(dp) function integral(nodes, elems, wtq, fq) result(r)
real(dp), intent(in) :: nodes(:, :)
integer, intent(in) :: elems(:, :)
real(dp), intent(in) :: wtq(:, :, :), fq(:, :, :, :)
real(dp) :: jacx, jacy, jacz, l(3), jac_det
integer :: e, Ne
Ne = size(elems, 2)
r = 0
do e = 1, Ne
    ! l is the diagonal vector:
    l = nodes(:, elems(7, e)) - nodes(:, elems(1, e))
    ! Assume rectangular shape:
    jacx = l(1)/2
    jacy = l(2)/2
    jacz = l(3)/2
    jac_det = abs(jacx*jacy*jacz)
    r = r + sum(fq(:, :, :, e) * jac_det * wtq)
end do
end function

end module

! ------------------------------------------------------------------------

module poisson3d_code

use types, only: dp
use feutils, only: phih, dphih
use fe_mesh, only: cartesian_mesh_3d, define_connect_tensor_3d, &
    c2fullc_3d, fe2quad_3d
use poisson3d_assembly, only: assemble_3d, sol_error, integral
use feutils, only: get_parent_nodes, get_parent_quad_pts_wts
use linalg, only: solve
use constants, only: pi
implicit none

contains

real(dp) function solve_poisson(Nex, Ney, Nez, p) result(error)
integer, intent(in) :: p

integer :: Nn, Ne
! nodes(:, i) are the (x,y) coordinates of the i-th mesh node
real(dp), allocatable :: nodes(:, :)
integer, allocatable :: elems(:, :) ! elems(:, i) are nodes of the i-th element
integer :: Nq, Nb
real(dp), allocatable :: xin(:), xiq(:), wtq(:), A(:, :), rhs(:), sol(:), &
        fullsol(:), solq(:, :, :, :), wtq3(:, :, :), phihq(:, :), dphihq(:, :)
integer, allocatable :: in(:, :, :, :), ib(:, :, :, :)
integer :: i, j, k
integer, intent(in) :: Nex, Ney, Nez

call cartesian_mesh_3d(Nex, Ney, Nez, &
    [0._dp, 0._dp, 0._dp], [1._dp, 1._dp, 1._dp], nodes, elems)
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
call define_connect_tensor_3d(Nex, Ney, Nez, p, 2, ib)
Nb = maxval(ib)
print *, "DOFs =", Nb
allocate(A(Nb, Nb), rhs(Nb), sol(Nb), fullsol(maxval(in)), solq(Nq, Nq, Nq, Ne))

call assemble_3d(xin, nodes, elems, ib, xiq, wtq3, phihq, dphihq, A, rhs)
sol = solve(A, rhs)
call c2fullc_3d(in, ib, sol, fullsol)
call fe2quad_3d(elems, xin, xiq, phihq, in, fullsol, solq)
error = sol_error(nodes, elems, xiq, wtq3, solq)
print *, "INTEGRAL:", integral(nodes, elems, wtq3, solq*solq)

end function

end module


! ------------------------------------------------------------------------


program poisson3d
use types, only: dp
use poisson3d_code, only: solve_poisson

real(dp) :: error

error = solve_poisson(1, 1, 1, 8)
print *, "L2 error:", error
end program
