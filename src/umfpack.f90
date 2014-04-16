module umfpack
use iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr
use types, only: dp
use utils, only: stop_error
implicit none

integer(c_int), parameter :: UMFPACK_OK = 0

integer(c_int), parameter :: UMFPACK_A = 0
integer(c_int), parameter :: UMFPACK_At = 1
integer(c_int), parameter :: UMFPACK_Aat = 2

integer(c_int), parameter :: UMFPACK_Pt_L = 3
integer(c_int), parameter :: UMFPACK_L = 4
integer(c_int), parameter :: UMFPACK_Lt_P = 5
integer(c_int), parameter :: UMFPACK_Lat_P = 6
integer(c_int), parameter :: UMFPACK_Lt = 7
integer(c_int), parameter :: UMFPACK_Lat = 8

integer(c_int), parameter :: UMFPACK_U_Qt = 9
integer(c_int), parameter :: UMFPACK_U = 10
integer(c_int), parameter :: UMFPACK_Q_Ut = 11
integer(c_int), parameter :: UMFPACK_Q_Uat = 12
integer(c_int), parameter :: UMFPACK_Ut = 13
integer(c_int), parameter :: UMFPACK_Uat = 14

interface

    integer(c_int) function umfpack_di_symbolic(m, n, Ap, Ai, Ax, Symbolic, &
        Control, Info) bind(c)
    import :: c_int, c_double, c_ptr
    integer(c_int), intent(in), value :: m, n
    integer(c_int), intent(in) :: Ap(n+1), Ai(*)
    real(c_double), intent(in) :: Ax(*)
    type(c_ptr), intent(out) :: Symbolic
    type(c_ptr), intent(in), value :: Control, Info
    end function

    integer(c_int) function umfpack_di_numeric(Ap, Ai, Ax, Symbolic, &
        Numeric, Control, Info) bind(c)
    import :: c_int, c_double, c_ptr
    integer(c_int), intent(in) :: Ap(*), Ai(*)
    real(c_double), intent(in) :: Ax(*)
    type(c_ptr), intent(in), value :: Symbolic
    type(c_ptr), intent(out) :: Numeric
    type(c_ptr), intent(in), value :: Control, Info
    end function

    integer(c_int) function umfpack_di_solve(sys, Ap, Ai, Ax, X, B, Numeric, &
        Control, Info) bind(c)
    import :: c_int, c_double, c_ptr
    integer(c_int), intent(in), value :: sys
    integer(c_int), intent(in) :: Ap(*), Ai(*)
    real(c_double), intent(in) :: Ax(*), B(*)
    real(c_double), intent(out) :: X(*)
    type(c_ptr), intent(in), value :: Numeric
    type(c_ptr), intent(in), value :: Control, Info
    end function

    subroutine umfpack_di_free_symbolic(Symbolic) bind(c)
    import :: c_ptr
    type(c_ptr), intent(inout) :: Symbolic
    end subroutine

    subroutine umfpack_di_free_numeric(Numeric) bind(c)
    import :: c_ptr
    type(c_ptr), intent(inout) :: Numeric
    end subroutine

end interface

type umfpack_numeric
    type(c_ptr) :: numeric
end type

contains

subroutine factorize(n, Ap, Ai, Ax, umfpack_data)
integer, intent(in) :: n, Ap(:), Ai(:)
real(dp), intent(in) :: Ax(:)
type(umfpack_numeric), intent(out) :: umfpack_data
type(c_ptr) :: Symbolic
integer :: s
s = umfpack_di_symbolic(n, n, Ap-1, Ai-1, Ax, Symbolic, c_null_ptr, c_null_ptr)
if (s /= UMFPACK_OK) call stop_error("umfpack_di_symbolic() failed")
s = umfpack_di_numeric(Ap-1, Ai-1, Ax, Symbolic, umfpack_data%numeric, &
    c_null_ptr, c_null_ptr);
if (s /= UMFPACK_OK) call stop_error("umfpack_di_numeric() failed")
call umfpack_di_free_symbolic(Symbolic);
end subroutine

subroutine solve(Ap, Ai, Ax, x, b, umfpack_data)
integer, intent(in) :: Ap(:), Ai(:)
real(dp), intent(in) :: Ax(:), b(:)
real(dp), intent(out) :: x(:)
type(umfpack_numeric), intent(in) :: umfpack_data
integer :: s
s = umfpack_di_solve(UMFPACK_A, Ap-1, Ai-1, Ax, x, b, umfpack_data%numeric, &
    c_null_ptr, c_null_ptr);
if (s /= UMFPACK_OK) call stop_error("umfpack_di_symbolic() failed")
end subroutine

subroutine free_data(umfpack_data)
type(umfpack_numeric), intent(inout) :: umfpack_data
call umfpack_di_free_numeric(umfpack_data%numeric)
end subroutine

end module
