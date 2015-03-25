module fe_mesh

! The quad is enumerated as follows:
!
! 4       3
!  +-----+
!  |     |
!  |     |
!  +-----+
! 1       2
!
! The hexahedron is enumerated as follows
!
! The bottom face:             The top face:
! 4       3                      8       7
!  +-----+                        +-----+
!  |     |                        |     |
!  |     |                        |     |
!  +-----+                        +-----+
! 1       2                      5       6
!

use types
use utils, only: stop_error, newunit, zeros
use quadrature, only: gauss_pts, gauss_wts, lobatto_wts, lobatto_pts
use solvers, only: solve_sym
use feutils, only: phih
use isolve, only: solve_cg
use utils, only: assert
implicit none
private
public cartesian_mesh_2d, cartesian_mesh_3d, define_connect_tensor_2d, &
    define_connect_tensor_3d, c2fullc_2d, c2fullc_3d, fe2quad_2d, fe2quad_3d, &
    vtk_save, get_element_3d, fe_eval_xyz, line_save, fe2quad_3d_lobatto, &
    quad2fe_3d, quad2fe_3d_lobatto, cartesian_mesh_3d_mask

contains

subroutine cartesian_mesh_2d(nx, ny, p1, p2, nodes, elems)
! Returns a uniform 2D cartesian mesh of nx x ny elements
!
! The lower left corner is p1(2), the upper right corner is p2(2).
integer, intent(in) :: nx, ny ! the number of elements in x and y-direction
real(dp), intent(in) :: p1(:), p2(:)
real(dp), allocatable, intent(out) :: nodes(:, :)
integer, allocatable, intent(out) :: elems(:, :)
integer :: Ne, Nn, i, j, idx
Ne = nx*ny
Nn = (nx+1)*(ny+1)
allocate(nodes(2, Nn), elems(4, Ne))
! Loop over nodes:
idx = 1
do i = 1, nx + 1
    do j = 1, ny + 1
        nodes(:, idx) = (p2-p1) * [i-1, j-1] / [nx, ny] + p1
        idx = idx + 1
    end do
end do
! Loop over elements:
idx = 1
elems = 0
do i = 1, nx
    do j = 1, ny
        elems(1, idx) = p(i  , j  )
        elems(2, idx) = p(i+1, j  )
        elems(3, idx) = p(i+1, j+1)
        elems(4, idx) = p(i  , j+1)
        idx = idx + 1
    end do
end do

contains

    integer pure function p(i, j)
    integer, intent(in) :: i, j
    p = (i-1)*(ny+1) + j
    end function

end subroutine

subroutine cartesian_mesh_3d(nx, ny, nz, p1, p2, nodes, elems)
! Returns a uniform 3D cartesian mesh of nx x ny x nz elements
!
! The lower left front corner is p1(3), the upper right back corner is p2(3).
integer, intent(in) :: nx, ny, nz ! number of elements in x, y and z-directions
real(dp), intent(in) :: p1(:), p2(:)
real(dp), allocatable, intent(out) :: nodes(:, :)
integer, allocatable, intent(out) :: elems(:, :)
integer :: Ne, Nn, i, j, k, idx
Ne = nx*ny*nz
Nn = (nx+1)*(ny+1)*(nz+1)
allocate(nodes(3, Nn), elems(8, Ne))
! Loop over nodes:
idx = 1
do i = 1, nx + 1
    do j = 1, ny + 1
        do k = 1, nz + 1
            nodes(:, idx) = (p2-p1) * [i-1, j-1, k-1] / [nx, ny, nz] + p1
            idx = idx + 1
        end do
    end do
end do
! Loop over elements:
idx = 1
elems = 0
do i = 1, nx
    do j = 1, ny
        do k = 1, nz
            elems(1, idx) = p(i  , j  , k  )
            elems(2, idx) = p(i+1, j  , k  )
            elems(3, idx) = p(i+1, j+1, k  )
            elems(4, idx) = p(i  , j+1, k  )
            elems(5, idx) = p(i  , j  , k+1)
            elems(6, idx) = p(i+1, j  , k+1)
            elems(7, idx) = p(i+1, j+1, k+1)
            elems(8, idx) = p(i  , j+1, k+1)
            idx = idx + 1
        end do
    end do
end do

contains

    integer pure function p(i, j, k)
    integer, intent(in) :: i, j, k
    p = (i-1)*(ny+1)*(nz+1) + (j-1)*(nz+1) + k
    end function

end subroutine


subroutine define_connect_tensor_2d(nex, ney, p, ibc, gn)
! 2D connectivity table for tensor-product order p elements
integer, intent(in) :: nex, ney ! Number of elements in x and y directions
integer, intent(in) :: p ! Polynomial order of elements
! Boundary condition: 1 = Neumann, 2 = Dirichlet, 3 = periodic
integer, intent(in) :: ibc
! gn(i, j, e) is the global node number of the local (i,j) node in e-th element
integer, allocatable, intent(out) :: gn(:, :, :)
integer :: nodes(nex*p+1, ney*p+1)
integer :: inode, ix, iy, iel, iex, iey
! Construct array of global nodes
! 0 = no associated basis function as in e.g. Dirichlet BCs
nodes = 0
inode = 0
do iy = 1, ney*p+1
    do ix = 1, nex*p+1
        if (ibc >= 2 .and. (ix == 1 .or. ix == nex*p+1 .or. &
                            iy == 1 .or. iy == ney*p+1)) cycle
        inode = inode + 1
        nodes(ix, iy) = inode
    end do
end do
if (ibc == 3) then
    ! Now we need to connect the basis functions on the opposite sites of the
    ! boundary.
    ! top-bottom sides
    do ix = 2, nex*p
        inode = inode + 1
        nodes(ix, 1)       = inode
        nodes(ix, ney*p+1) = inode
    end do
    ! left-right sides
    do iy = 2, ney*p
        inode = inode + 1
        nodes(1, iy)       = inode
        nodes(nex*p+1, iy) = inode
    end do
    ! Corners
    inode = inode + 1
    nodes(1, 1)             = inode
    nodes(1, ney*p+1)       = inode
    nodes(nex*p+1, 1)       = inode
    nodes(nex*p+1, ney*p+1) = inode
end if
! Construct connectivity table of global nodes in each element
allocate(gn(p+1, p+1, nex*ney))
iel = 0
do iex = 1, nex
    do iey = 1, ney
        iel = iel + 1 ! Element number
        ix = (iex-1)*p+1 ! Lower left corner
        iy = (iey-1)*p+1
        ! Get global node numbers in element
        gn(:, :, iel) = nodes(ix:ix+p, iy:iy+p)
    end do
end do
end subroutine

subroutine define_connect_tensor_3d(nex, ney, nez, p, ibc, gn)
! 3D connectivity table for tensor-product order p elements
integer, intent(in) :: nex, ney, nez ! Number of elements in x, y, z directions
integer, intent(in) :: p ! Polynomial order of elements
! Boundary condition: 1 = Neumann, 2 = Dirichlet, 3 = periodic
integer, intent(in) :: ibc
! gn(i, j, k, e) is the global node number of the local (i,j,k) node
! in the e-th element:
integer, allocatable, intent(out) :: gn(:, :, :, :)
integer :: nodes(nex*p+1, ney*p+1, nez*p+1)
integer :: inode, ix, iy, iz, iel, iex, iey, iez
! Construct array of global nodes
! 0 = no associated basis function as in e.g. Dirichlet BCs
nodes = 0
inode = 0
do iz = 1, nez*p+1
    do iy = 1, ney*p+1
        do ix = 1, nex*p+1
            if (ibc >= 2 .and. (ix == 1 .or. ix == nex*p+1 .or. &
                iy == 1 .or. iy == ney*p+1 .or. &
                iz == 1 .or. iz == nez*p+1)) cycle
            inode = inode + 1
            nodes(ix, iy, iz) = inode
        end do
    end do
end do
if (ibc == 3) then
    ! Now we need to connect the basis functions on the opposite sites of the
    ! boundary.
    ! top-bottom faces (middle part)
    do ix = 2, nex*p
        do iy = 2, ney*p
            inode = inode + 1
            nodes(ix, iy, 1)       = inode
            nodes(ix, iy, nez*p+1) = inode
        end do
    end do
    ! top-bottom faces (4 edges connecting them)
    do iz = 2, nez*p
        inode = inode + 1
        nodes(1, 1, iz)             = inode
        nodes(1, ney*p+1, iz)       = inode
        nodes(nex*p+1, 1, iz)       = inode
        nodes(nex*p+1, ney*p+1, iz) = inode
    end do
    ! left-right faces (middle part)
    do iy = 2, ney*p
        do iz = 2, nez*p
            inode = inode + 1
            nodes(1, iy, iz)       = inode
            nodes(nex*p+1, iy, iz) = inode
        end do
    end do
    ! left-right faces (4 edges connecting them)
    do ix = 2, nex*p
        inode = inode + 1
        nodes(ix, 1, 1)             = inode
        nodes(ix, 1, nez*p+1)       = inode
        nodes(ix, ney*p+1, 1)       = inode
        nodes(ix, ney*p+1, nez*p+1) = inode
    end do
    ! front-back faces (middle part)
    do ix = 2, nex*p
        do iz = 2, nez*p
            inode = inode + 1
            nodes(ix, 1, iz)       = inode
            nodes(ix, ney*p+1, iz) = inode
        end do
    end do
    ! front-back faces (4 edges connecting them)
    do iy = 2, ney*p
        inode = inode + 1
        nodes(1, iy, 1)             = inode
        nodes(1, iy, nez*p+1)       = inode
        nodes(nex*p+1, iy, 1)       = inode
        nodes(nex*p+1, iy, nez*p+1) = inode
    end do
    ! Corners (vertices)
    inode = inode + 1
    nodes(1, 1, 1)                   = inode
    nodes(1, ney*p+1, 1)             = inode
    nodes(nex*p+1, 1, 1)             = inode
    nodes(nex*p+1, ney*p+1, 1)       = inode
    nodes(1, 1, nez*p+1)             = inode
    nodes(1, ney*p+1, nez*p+1)       = inode
    nodes(nex*p+1, 1, nez*p+1)       = inode
    nodes(nex*p+1, ney*p+1, nez*p+1) = inode
end if
! Construct connectivity table of global nodes in each element
allocate(gn(p+1, p+1, p+1, Nex*Ney*Nez))
iel = 0
do iex = 1, Nex
    do iey = 1, Ney
        do iez = 1, Nez
            iel = iel + 1 ! Element number
            ix = (iex-1)*p+1 ! Lower left front corner
            iy = (iey-1)*p+1
            iz = (iez-1)*p+1
            ! Get global node numbers in element
            gn(:, :, :, iel) = nodes(ix:ix+p, iy:iy+p, iz:iz+p)
        end do
    end do
end do
end subroutine


subroutine c2fullc_2d(in, ib, c, fullc)
! Converts FE coefficient vector to full coefficient vector
! It puts 0 for Dirichlet boundary conditions (ib==0), otherwise it just copies
! the coefficients.
integer, intent(in) :: in(:, :, :)
integer, intent(in) :: ib(:, :, :)
real(dp), intent(in) :: c(:) ! coefficient vector with regards to ib
real(dp), intent(out) :: fullc(:) ! full coefficients vector with regards to in
integer :: e, i, j
do e = 1, size(in, 3)
    do i = 1, size(in, 1)
        do j = 1, size(in, 2)
            if (ib(i, j, e) == 0) then
                fullc(in(i, j, e)) = 0 ! Dirichlet
            else
                fullc(in(i, j, e)) = c(ib(i, j, e))
            end if
        end do
    end do
end do
end subroutine

subroutine c2fullc_3d(in, ib, c, fullc)
! Converts FE coefficient vector to full coefficient vector
! It puts 0 for Dirichlet boundary conditions (ib==0), otherwise it just copies
! the coefficients.
integer, intent(in) :: in(:, :, :, :)
integer, intent(in) :: ib(:, :, :, :)
real(dp), intent(in) :: c(:) ! coefficient vector with regards to ib
real(dp), intent(out) :: fullc(:) ! full coefficients vector with regards to in
integer :: e, i, j, k
do e = 1, size(in, 4)
    do i = 1, size(in, 1)
        do j = 1, size(in, 2)
            do k = 1, size(in, 3)
                if (ib(i, j, k, e) == 0) then
                    fullc(in(i, j, k, e)) = 0 ! Dirichlet
                else
                    fullc(in(i, j, k, e)) = c(ib(i, j, k, e))
                end if
            end do
        end do
    end do
end do
end subroutine


subroutine fe2quad_2d(elems, xin, xiq, phihq, in, fullu, uq)
! Transforms fullu from FE-coefficient to quadrature-grid representation.
! fullu is a full FE coefficient vector, having values for all nodes in the
! mesh, including domain-boundary nodes.
integer, intent(in) :: elems(:, :)
real(dp), intent(in) :: xin(:)
real(dp), intent(in) :: xiq(:)
real(dp), intent(in) :: phihq(:, :)
integer, intent(in) :: in(:, :, :)
real(dp), intent(in) :: fullu(:)
real(dp), intent(out) :: uq(:, :, :)
integer :: ie, ilnx, ilny, iqx, iqy

! evaluate at quad points in each element
do ie = 1, size(elems, 2)
    uq(:, :, ie) = 0
    do ilnx = 1, size(xin)
    do ilny = 1, size(xin)
        do iqx = 1, size(xiq)
        do iqy = 1, size(xiq)
            uq(iqx, iqy, ie) = uq(iqx, iqy, ie) + &
                fullu(in(ilnx, ilny, ie)) * phihq(iqx, ilnx) * phihq(iqy, ilny)
        end do
        end do
    end do
    end do
end do
end subroutine

subroutine fe2quad_3d(elems, xin, xiq, phihq, in, fullu, uq)
! Transforms fullu from FE-coefficient to quadrature-grid representation.
! fullu is a full FE coefficient vector, having values for all nodes in the
! mesh, including domain-boundary nodes.
integer, intent(in) :: elems(:, :)
real(dp), intent(in) :: xin(:)
real(dp), intent(in) :: xiq(:)
real(dp), intent(in) :: phihq(:, :)
integer, intent(in) :: in(:, :, :, :)
real(dp), intent(in) :: fullu(:)
real(dp), intent(out) :: uq(:, :, :, :)
integer :: ie, ilnx, ilny, ilnz, iqx, iqy, iqz

! evaluate at quad points in each element
do ie = 1, size(elems, 2)
    uq(:, :, :, ie) = 0
    do ilnx = 1, size(xin)
    do ilny = 1, size(xin)
    do ilnz = 1, size(xin)
        do iqx = 1, size(xiq)
        do iqy = 1, size(xiq)
        do iqz = 1, size(xiq)
            uq(iqx, iqy, iqz, ie) = uq(iqx, iqy, iqz, ie) + &
                fullu(in(ilnx, ilny, ilnz, ie)) * &
                phihq(iqx, ilnx) * phihq(iqy, ilny) * phihq(iqz, ilnz)
        end do
        end do
        end do
    end do
    end do
    end do
end do
end subroutine

subroutine quad2fe_3d(Ne, Nb, p, jac_det, wtq3, Sp, Sj, Sx, phi_v, in, ib, &
        uq, fullu)
! Transforms quadrature-grid representation to FE-coefficient (fullu).
! fullu is a full FE coefficient vector, having values for all nodes in the
! mesh, including domain-boundary nodes.
! Assumes the same 'jac_det' for each finite element.
use poisson3d_assembly, only: assemble_3d_coo_rhs
integer, intent(in) :: Ne, Nb, p
real(dp), intent(in) :: jac_det ! |J| for each element (must be constant)
real(dp), intent(in) :: wtq3(:, :, :)
integer, intent(in) :: Sp(:), Sj(:)
real(dp), intent(in) :: Sx(:)
real(dp), intent(in) :: phi_v(:, :, :, :, :, :)
real(dp), intent(in) :: uq(:, :, :, :)
integer, intent(in) :: in(:, :, :, :), ib(:, :, :, :)
real(dp), intent(out) :: fullu(:)
real(dp) :: rhs(Nb), sol(Nb)
call assemble_3d_coo_rhs(Ne, p, uq, jac_det, wtq3, ib, phi_v, rhs)
sol = solve_cg(Sp, Sj, Sx, rhs, zeros(size(rhs)), 1e-12_dp, 400)
call c2fullc_3d(in, ib, sol, fullu)
end subroutine

subroutine quad2fe_3d_lobatto(Ne, p, in, uq, fullu)
! Transforms quadrature-grid representation to FE-coefficient (fullu).
! fullu is a full FE coefficient vector, having values for all nodes in the
! mesh, including domain-boundary nodes.
! Assumes the same 'jac_det' for each finite element and spectral elements.
integer, intent(in) :: Ne, p
real(dp), intent(in) :: uq(:, :, :, :)
integer, intent(in) :: in(:, :, :, :)
real(dp), intent(out) :: fullu(:)
integer :: ie, iqx, iqy, iqz
! The quad points are the coefficients for spectral elements
do ie = 1, Ne
    do iqz = 1, p+1
    do iqy = 1, p+1
    do iqx = 1, p+1
        fullu(in(iqx, iqy, iqz, ie)) = uq(iqx, iqy, iqz, ie)
    end do
    end do
    end do
end do
end subroutine

subroutine fe2quad_3d_lobatto(elems, xiq, in, fullu, uq)
integer, intent(in) :: elems(:, :)
real(dp), intent(in) :: xiq(:)
integer, intent(in) :: in(:, :, :, :)
real(dp), intent(in) :: fullu(:)
real(dp), intent(out) :: uq(:, :, :, :)
integer :: ie, iqx, iqy, iqz
! evaluate at quad points in each element
do ie = 1, size(elems, 2)
    do iqz = 1, size(xiq)
    do iqy = 1, size(xiq)
    do iqx = 1, size(xiq)
        uq(iqx, iqy, iqz, ie) = fullu(in(iqx, iqy, iqz, ie))
    end do
    end do
    end do
end do
end subroutine

subroutine vtk_save(filename, nx, ny, nz, nodes, elems, xiq, fq)
character(len=*), intent(in) :: filename
integer, intent(in) :: nx, ny, nz
real(dp), intent(in):: nodes(:, :)
integer, intent(in) :: elems(:, :)
real(dp), intent(in) :: xiq(:)
real(dp), intent(in) :: fq(:, :, :, :)
integer :: Ne, e, iqx, iqy, iqz, i, j, k, Nq
real(dp), dimension(size(xiq)) :: x, y, z, xp, yp, zp
real(dp) :: qdata(nx*size(xiq), ny*size(xiq), nz*size(xiq))
real(dp) :: lx, ly, lz
real(dp) :: jacx, jacy, jacz, jac_det
integer :: u

Nq = size(xiq)

open(newunit(u), file=filename, status="replace")
write(u, "(a)") "# vtk DataFile Version 2.0"
write(u, "(a)") "Function at quadrature points"
write(u, "(a)") "ASCII"
write(u, "(a)") "DATASET RECTILINEAR_GRID"
write(u, "(a, i8, i8, i8)") "DIMENSIONS", nx*Nq, ny*Nq, nz*Nq

Ne = size(elems, 2)
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
! Note: this depends on the exact order 'elems' is constructed
write(u, "(a, ' ', i8, ' ', a)") "X_COORDINATES", nx*Nq, "float"
do e = 1, nx*ny*nz, ny*nz
    x = xp + nodes(1, elems(1, e))
    do iqz = 1, Nq
        write(u, *) x(iqz)
    end do
end do
write(u, "(a, ' ', i8, ' ', a)") "Y_COORDINATES", ny*Nq, "float"
do e = 1, ny*nz, nz
    y = yp + nodes(2, elems(1, e))
    do iqz = 1, Nq
        write(u, *) y(iqz)
    end do
end do
write(u, "(a, ' ', i8, ' ', a)") "Z_COORDINATES", nz*Nq, "float"
do e = 1, nz
    z = zp + nodes(3, elems(1, e))
    do iqz = 1, Nq
        write(u, *) z(iqz)
    end do
end do

e = 1
do i = 1, nx
do j = 1, ny
do k = 1, nz
    do iqx = 1, Nq
    do iqy = 1, Nq
    do iqz = 1, Nq
        qdata((i-1)*Nq + iqx, (j-1)*Nq + iqy, (k-1)*Nq + iqz) = &
            fq(iqx, iqy, iqz, e)
    end do
    end do
    end do
    e = e + 1
end do
end do
end do

write(u, "(a, ' ', i10)") "POINT_DATA", size(qdata)
write(u, "(a)") "SCALARS MyFunction float"
write(u, "(a)") "LOOKUP_TABLE default"
write(u, *) qdata

close(u)
end subroutine

subroutine line_save(filename, xin, nodes, elems, in, fullu, x1, x2, n)
character(len=*), intent(in) :: filename
real(dp), intent(in):: xin(:), nodes(:, :), fullu(:), x1(:), x2(:)
integer, intent(in):: elems(:, :), in(:, :, :, :), n
integer :: i, u
real(dp) :: x(3), val
open(newunit(u), file=filename, status="replace")
do i = 1, n
    x = x1 + (i-1) * (x2-x1) / (n-1)
    val = fe_eval_xyz(xin, nodes, elems, in, fullu, x)
    write(u, *) x, val
end do
close(u)
end subroutine

real(dp) function fe_eval_xyz(xin, nodes, elems, in, fullu, x) result(val)
real(dp), intent(in):: xin(:), nodes(:, :), x(:), fullu(:)
integer, intent(in):: elems(:, :), in(:, :, :, :)
real(dp), dimension(3) :: l, jac, xi, xa, xb
integer :: nx, ny, nz, p, e

e = get_element_3d(nodes, elems, x)

! xa, xb are two opposite corners of the hexahedron
xa = nodes(:, elems(1, e))
xb = nodes(:, elems(7, e))
! l is the diagonal vector:
l = xb - xa
! Assume hexahedral shape:
jac = l / 2
! xi is in the reference domain [-1, 1]^3
xi = (x-xa) / jac - 1

p = size(xin) - 1
val = 0
do nz = 1, p+1
do ny = 1, p+1
do nx = 1, p+1
    val = val + fullu(in(nx, ny, nz, e)) * &
        phih(xin, nx, xi(1)) * phih(xin, ny, xi(2)) * phih(xin, nz, xi(3))
end do
end do
end do
end function

integer function get_element_3d(nodes, elems, x) result(e)
! Returns the element that containst the point 'x'
real(dp), intent(in) :: nodes(:, :), x(:)
integer, intent(in) :: elems(:, :)
real(dp) :: xa(3), xb(3)
integer :: i, Ne
Ne = size(elems, 2)
do i = 1, Ne
    ! xa, xb are two opposite corners of the hexahedron
    xa = nodes(:, elems(1, i))
    xb = nodes(:, elems(7, i))
    if (point_in_hex(xa, xb, x)) then
        e = i
        return
    end if
end do
e = 0
call stop_error("The point 'x' is not inside the mesh.")
end function

logical pure function point_in_hex(p1, p2, x) result(r)
! Returns .true. if the point 'x' is in the hexahedron specified by the two
! opposite corners p1 (lower, left, front) and p2 (upper, right, back).
real(dp), intent(in) :: p1(:), p2(:), x(:)
r = all(p1 <= x) .and. all(x <= p2)
end function


subroutine cartesian_mesh_3d_mask(myid, nx, ny, nz, nsubx, nsuby, nsubz, &
        mask_elems)
! Returns a mask for a subdomain "myid" within a uniform 3D cartesian mesh of
! (nx, ny, nz) elements, where subdomains have (nsubx, nsuby, nsubz) elements.
! Total number of elements in the domain in each direction
integer, intent(in) :: nx, ny, nz
! number of subdomains elements in x, y and z-directions
integer, intent(in) :: nsubx, nsuby, nsubz
! Processor id = 0, 1, 2, ...
integer, intent(in) :: myid
integer, allocatable, intent(out) :: mask_elems(:)
integer :: Ne, i, j, k, idx
integer :: nesubx, nesuby, nesubz
integer :: iminmax(6)
call assert(modulo(nx, nsubx) == 0)
call assert(modulo(ny, nsuby) == 0)
call assert(modulo(nz, nsubz) == 0)
Ne = nx*ny*nz
nesubx = nx / nsubx
nesuby = ny / nsuby
nesubz = nz / nsubz
allocate(mask_elems(Ne))

iminmax = get_minmax(myid+1, nsubx, nsuby, nsubz, nesubx, nesuby, nesubz)

mask_elems = 0
idx = 1
do i = 1, nx
    do j = 1, ny
        do k = 1, nz
            if     (i >= iminmax(1) .and. i <= iminmax(2) .and. &
                    j >= iminmax(3) .and. j <= iminmax(4) .and. &
                    k >= iminmax(5) .and. k <= iminmax(6)) then
                mask_elems(idx) = 1
            end if
            idx = idx + 1
        end do
    end do
end do
end subroutine

function get_minmax(id, nsubx, nsuby, nsubz, nesubx, nesuby, nesubz) result(r)
! Calculates element indices belonging to the subdomain 'id'
integer, intent(in) :: id ! The subdomain ID = 1, 2, 3, ...
! number of subdomains in each direction:
integer, intent(in) :: nsubx, nsuby, nsubz
! number of elements in a subdomain in each direction:
integer, intent(in) :: nesubx, nesuby, nesubz
integer :: r(6) ! r = [iminx, imaxx, iminy, imaxy, iminz, imaxz]
integer :: i, j, k, idx
idx = 1
do i = 1, nsubx
    do j = 1, nsuby
        do k = 1, nsubz
            if (idx == id) then
                r = [(i-1)*nesubx + 1, i*nesubx, &
                     (j-1)*nesuby + 1, j*nesuby, &
                     (k-1)*nesubz + 1, k*nesubz]
                return
            end if
            idx = idx + 1
        end do
    end do
end do
call stop_error("get_minmax(): 'id' not found")
end function

end module
