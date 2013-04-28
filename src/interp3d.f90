module interp3d
use types, only: dp
implicit none
private
public trilinear

contains

real(dp) function trilinear(x, p1, p2, values) result(r)
! Returns trilinear interpolation of a point x0 using uniform data values
real(dp), intent(in) :: x(:) ! The 3D coordinates of a point to interpolate
real(dp), intent(in) :: p1(:) ! The lower left front corner is p1(3)
real(dp), intent(in) :: p2(:) ! The upper right back corner is p2(3)
real(dp), intent(in) :: values(:, :, :) ! Values on a uniform 3D grid
real(dp) :: x0(3)
integer :: ijk(3), nelem(3)
! Number of elements in each direction:
nelem = shape(values) - 1
! Transform to box [0,1] x [0,1] x [0,1] * nelem, that is, make the
! element length exactly 1 in each direction
x0 = nelem * (x - p1) / (p2 - p1)
ijk = int(x0)+1 ! indices of the nearest vertex
where (ijk > nelem) ijk = nelem
where (ijk < 1) ijk = 1
r = trilinear_unitbox(x0-ijk+1, values(ijk(1):ijk(1)+1, ijk(2):ijk(2)+1, &
        ijk(3):ijk(3)+1))
end function

real(dp) function trilinear_unitbox(xyz, V) result(Vxyz)
! Returns trilinear interpolation of a point (x, y, z) in a unit box
real(dp), intent(in) :: xyz(3)  ! The point in a unit box [0,1] x [0,1] x [0,1]
! V(i, j, k) ... value at vertex i, j, k of the box, where i,j,k = 1,2:
real(dp), intent(in) :: V(:, :, :)
real(dp) :: x, y, z
x = xyz(1)
y = xyz(2)
z = xyz(3)
Vxyz = &
    V(1, 1, 1) * (1-x)*(1-y)*(1-z) + &
    V(2, 1, 1) *   x  *(1-y)*(1-z) + &
    V(1, 2, 1) * (1-x)*  y  *(1-z) + &
    V(1, 1, 2) * (1-x)*(1-y)*  z   + &
    V(1, 2, 2) * (1-x)*  y  *  z   + &
    V(2, 1, 2) *   x  *(1-y)*  z   + &
    V(2, 2, 1) *   x  *  y  *(1-z) + &
    V(2, 2, 2) *   x  *  y  *  z
end function

end module
