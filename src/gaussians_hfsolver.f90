module gaussians

! Two particle Gaussian integrals are calculated using hfsolver

use types, only: dp
use qc, only: getInts
implicit none
private
public getints2

contains

subroutine getints2(nprim, istart, center, power, coef, alpha, int2)
integer, intent(in) :: nprim(:), istart(:)
real(dp), intent(in) :: center(:, :)
integer, intent(in) :: power(:, :)
real(dp), intent(in) :: coef(:), alpha(:)
real(dp), intent(out) :: int2(:)

integer, allocatable :: lpower(:), mpower(:), npower(:)
real(dp), allocatable :: xcenter(:), ycenter(:), zcenter(:)
integer :: n

n = size(nprim)
allocate(lpower(n), mpower(n), npower(n), xcenter(n), ycenter(n), zcenter(n))
lpower = power(1, :)
mpower = power(2, :)
npower = power(3, :)
xcenter = center(1, :)
ycenter = center(2, :)
zcenter = center(3, :)
call getInts(n, nprim, istart-1, xcenter, ycenter, zcenter, &
        lpower, mpower, npower, size(coef), coef, alpha, int2)
end subroutine

end module
