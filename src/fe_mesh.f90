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
use types
use utils, only: stop_error
use quadrature, only: gauss_pts, gauss_wts, lobatto_wts, lobatto_pts
use solvers, only: solve_sym
implicit none
private
public cartesian_mesh_2d, cartesian_mesh_3d, define_connect_tensor_2d, &
    c2fullc_2d, fe2quad_2d

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
integer, intent(in) :: ibc ! Boundary condition: 1 = Neumann, 2 = Dirichlet
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
        if (ibc == 2 .and. (ix == 1 .or. ix == nex*p+1 .or. &
                            iy == 1 .or. iy == ney*p+1)) cycle
        inode = inode + 1
        nodes(ix, iy) = inode
    end do
end do
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
integer, intent(in) :: ibc ! Boundary condition: 1 = Neumann, 2 = Dirichlet
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
            if (ibc == 2 .and. (ix == 1 .or. ix == nex*p+1 .or. &
                iy == 1 .or. iy == ney*p+1 .or. &
                iz == 1 .or. iz == nez*p+1)) cycle
            inode = inode + 1
            nodes(ix, iy, iz) = inode
        end do
    end do
end do
! Construct connectivity table of global nodes in each element
allocate(gn(p+1, p+1, p+1, nex*ney))
iel = 0
do iex = 1, nex
    do iey = 1, ney
        do iez = 1, nez
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

end module
