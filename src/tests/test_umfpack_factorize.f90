program test_umfpack_factorize
use types, only: dp
use umfpack, only: factorize, solve, free_data, umfpack_numeric
use utils, only: assert
implicit none

integer, parameter :: n = 5
integer, parameter :: Ap(*) = [1, 3, 6, 10, 11, 13]
integer, parameter :: Ai(*) = [1, 2, 1, 3, 5, 2, 3, 4, 5, 3, 2, 5]
real(dp), parameter :: Ax(*) = [2._dp, 3._dp, 3._dp, -1._dp, 4._dp, 4._dp, &
        -3._dp, 1._dp, 2._dp, 2._dp, 6._dp, 1._dp]
real(dp), parameter :: b(*) = [8._dp, 45._dp, -3._dp, 3._dp, 19._dp]
real(dp), parameter :: b2(*) = [-8._dp, -45._dp, 3._dp, -3._dp, -19._dp]
real(dp) :: x(n)
integer :: i

type(umfpack_numeric) :: d
call factorize(n, Ap, Ai, Ax, d)
call solve(Ap, Ai, Ax, x, b, d)

print *, x
do i = 1, n
    print "('x[', i0, '] = ', f0.6)", i, x(i)
    call assert(abs(x(i) - i) < 1e-10_dp);
end do

call solve(Ap, Ai, Ax, x, b2, d)

print *, x
do i = 1, n
    print "('x[', i0, '] = ', f0.6)", i, x(i)
    call assert(abs(x(i) - (-i)) < 1e-10_dp);
end do

call free_data(d);
end program
