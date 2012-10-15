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
public cartesian_mesh_2d, cartesian_mesh_3d

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

end module
