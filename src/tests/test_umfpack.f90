program test_umfpack
use iso_c_binding, only: c_ptr, c_null_ptr
use types, only: dp
use umfpack, only: umfpack_di_symbolic, umfpack_di_numeric, umfpack_di_solve, &
        umfpack_di_free_symbolic, umfpack_di_free_numeric, UMFPACK_A
use utils, only: assert
implicit none

integer, parameter :: n = 5
integer, parameter :: Ap(*) = [0, 2, 5, 9, 10, 12]
integer, parameter :: Ai(*) = [0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4]
real(dp), parameter :: Ax(*) = [2._dp, 3._dp, 3._dp, -1._dp, 4._dp, 4._dp, &
        -3._dp, 1._dp, 2._dp, 2._dp, 6._dp, 1._dp]
real(dp), parameter :: b(*) = [8._dp, 45._dp, -3._dp, 3._dp, 19._dp]
real(dp) :: x(n)
integer :: i, s

type(c_ptr) :: Symbolic, Numeric
s = umfpack_di_symbolic(n, n, Ap, Ai, Ax, Symbolic, c_null_ptr, c_null_ptr)
s = umfpack_di_numeric(Ap, Ai, Ax, Symbolic, Numeric, c_null_ptr, c_null_ptr)
call umfpack_di_free_symbolic(Symbolic)
s = umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, c_null_ptr, &
    c_null_ptr)
call umfpack_di_free_numeric(Numeric)

print *, x
do i = 1, n
    print "('x[', i0, '] = ', f0.6)", i, x(i)
    call assert(abs(x(i) - i) < 1e-10_dp)
end do
end program
