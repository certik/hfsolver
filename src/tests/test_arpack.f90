program test_arpack
use arpack, only: eig
use types, only: dp
use utils, only: stop_error
implicit none

integer :: n, nev, ncv, i
real(dp), allocatable :: d(:), v(:,:)
real(dp), parameter :: e_ref(4) =  [4.8308300260037749_dp, &
    5.3097214678905722_dp, 5.6825070656623593_dp, 5.9189859472289967_dp]


n = 10
nev = 4
ncv = 10
allocate(v(n,ncv), d(ncv))
call eig(n, nev, ncv, "LM", av, d, v)
do i = 1, nev
    print *, d(i), abs(d(i) - e_ref(i))
    if (abs(d(i) - e_ref(i)) > 1e-14_dp) call stop_error("Error in ref. check")
end do

contains

  subroutine av(x, y)
  ! Compute y = A*x
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: y(:)
  real(dp) :: dd, dl, du
  integer :: j
  dd = 4
  dl = -1
  du = -1
  y(1) = dd*x(1) + du*x(2)
  do j = 2, size(y)-1
     y(j) = dl*x(j-1) + dd*x(j) + du*x(j+1)
  end do
  y(size(y)) = dl*x(size(y)-1) + dd*x(size(y))
  end

end
