module schroed_assembly
use types, only: dp
use sorting, only: sort
implicit none
private
public assemble_3d

contains

real(dp) elemental function f(x, y, z)
real(dp), intent(in) :: x, y, z
real(dp) :: V0, r0, r, Xion(3, 2)
integer :: n
V0 = 16
r0 = 0.5_dp
f = 0
Xion(:,1) = 5._dp/2
Xion(:,2) = 5._dp/2
Xion(1,:) = Xion(1,:) + [-1, 1]
do n = 1, 2
    r = sqrt(sum(([x,y,z]-Xion(:,n))**2))
    f = f - V0 * exp(-r**2/r0**2)
end do
end function

subroutine assemble_3d(xin, nodes, elems, ib, xiq, wtq, phihq, dphihq, Am, Bm)
! Assemble on a 2D rectangular uniform mesh
real(dp), intent(in):: xin(:), nodes(:, :), xiq(:), wtq(:, :, :), &
    phihq(:, :), dphihq(:, :)
integer, intent(in):: elems(:, :), ib(:, :, :, :)
real(dp), intent(out):: Am(:,:), Bm(:, :)
real(dp), dimension(size(xiq), size(xiq), size(xiq), &
    size(xin), size(xin), size(xin)) :: phi_v, phi_dx, phi_dy, phi_dz
real(dp), dimension(size(xiq), size(xiq), size(xiq), &
    size(xiq), size(xiq), size(xiq)) :: Am_loc, Bm_loc
real(dp), dimension(size(xiq), size(xiq), size(xiq)) :: fq
real(dp), dimension(size(xiq)) :: x, y, z, xp, yp, zp
integer :: Ne, p, e, i, j, iqx, iqy, iqz
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
Am=0; Bm=0
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
! Precalculate element matrices:
do bz = 1, p+1
do by = 1, p+1
do bx = 1, p+1
    do az = 1, p+1
    do ay = 1, p+1
    do ax = 1, p+1
        Am_loc(ax, ay, az, bx, by, bz) = sum(( &
            phi_dx(:, :, :, ax, ay, az)*phi_dx(:, :, :, bx, by, bz) + &
            phi_dy(:, :, :, ax, ay, az)*phi_dy(:, :, :, bx, by, bz) + &
            phi_dz(:, :, :, ax, ay, az)*phi_dz(:, :, :, bx, by, bz)) &
            * jac_det * wtq) / 2
        Bm_loc(ax, ay, az, bx, by, bz) = sum(( &
            phi_v(:, :, :, ax, ay, az) * phi_v(:, :, :, bx, by, bz) &
            * jac_det * wtq))
    end do
    end do
    end do
end do
end do
end do
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
            Am(i,j) = Am(i,j) + Am_loc(ax, ay, az, bx, by, bz)
            Am(i,j) = Am(i,j) + sum(fq * &
                phi_v(:, :, :, ax, ay, az)*phi_v(:, :, :, bx, by, bz) &
                * jac_det * wtq)
            Bm(i,j) = Bm(i,j) + Bm_loc(ax, ay, az, bx, by, bz)
        end do
        end do
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


!------------------------------------------------------------------------------

program schroed3d

use types, only: dp
use fe_mesh, only: cartesian_mesh_3d, define_connect_tensor_3d
use schroed_assembly, only: assemble_3d
use feutils, only: get_parent_nodes, get_parent_quad_pts_wts, phih, dphih
use linalg, only: eigh
use linalg_feast, only: eigh_feast => eigh
use arpack, only: eig
implicit none

integer :: Nn, Ne
! nodes(:, i) are the (x,y) coordinates of the i-th mesh node
real(dp), allocatable :: nodes(:, :)
integer, allocatable :: elems(:, :) ! elems(:, i) are nodes of the i-th element
integer :: Nq, p, Nb
real(dp), allocatable :: xin(:), xiq(:), wtq(:), A(:, :), B(:, :), c(:, :), &
    lam(:), wtq3(:, :, :), phihq(:, :), dphihq(:, :), v(:,:), d(:)
integer, allocatable :: ib(:, :, :, :), in(:, :, :, :)
integer :: i, j, k, Nex, Ney, Nez, Neig, solver_type, M0, nev, ncv

Nex = 4
Ney = 4
Nez = 4
p = 2
Nq = p+1

call cartesian_mesh_3d(Nex, Ney, Nez, &
    [0._dp, 0._dp, 0._dp], [6._dp, 5._dp, 5._dp], &
    nodes, elems)
Nn = size(nodes, 2)
Ne = size(elems, 2)

print *, "Number of nodes:", Nn
print *, "Number of elements:", Ne

allocate(xin(p+1))
call get_parent_nodes(2, p, xin)
allocate(xiq(Nq), wtq(Nq), wtq3(Nq, Nq, Nq))
call get_parent_quad_pts_wts(2, Nq, xiq, wtq)
forall(i=1:Nq, j=1:Nq, k=1:Nq) wtq3(i, j, k) = wtq(i)*wtq(j)*wtq(k)
allocate(phihq(size(xiq), size(xin)))
allocate(dphihq(size(xiq), size(xin)))
! Tabulate parent basis at quadrature points
forall(i=1:size(xiq), j=1:size(xin))  phihq(i, j) =  phih(xin, j, xiq(i))
forall(i=1:size(xiq), j=1:size(xin)) dphihq(i, j) = dphih(xin, j, xiq(i))

call define_connect_tensor_3d(Nex, Ney, Nez, p, 1, in)
call define_connect_tensor_3d(Nex, Ney, Nez, p, 3, ib)
Nb = maxval(ib)
print *, "p =", p
print *, "DOFs =", Nb
allocate(A(Nb, Nb), B(Nb, Nb))

print *, "Assembling..."
call assemble_3d(xin, nodes, elems, ib, xiq, wtq3, phihq, dphihq, A, B)
!print *, "Solving..."
!solver_type = 1
!select case(solver_type)
!    case (1)
Neig = Nb
allocate(c(Nb, Neig), lam(Neig))
!        call eigh(A, B, lam, c)
!    case (2)
!        M0 = 15
!        call eigh_feast(A, B, 0._dp, 5._dp, M0, lam, c)
!        Neig = size(lam)
!end select
!print *, "Eigenvalues:"
Neig = 10
!do i = 1, Neig
!    print *, i, lam(i)
!end do
nev = 1
ncv = 20
allocate(v(Nb,ncv), d(ncv))
call eig(Nb, nev, ncv, "SA", av, d, v)
do i = 1, nev
    print *, i, d(i)
end do

do i = 1, Nb
    B(i,i) = 1 / sqrt(B(i,i))
end do
print *, "Mul"
A = matmul(matmul(B, A), B)
print *, "Solve"
call eigh(A, lam, c)
do i = 1, Neig
    print *, i, lam(i)
end do

contains

    subroutine av(x, y)
    ! Compute y = A*x
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    do i = 1, Nb
        y(i) = x(i) / sqrt(B(i,i))
    end do
    y = matmul(A, y)
    do i = 1, Nb
        y(i) = y(i) / sqrt(B(i,i))
    end do
    end
end program
