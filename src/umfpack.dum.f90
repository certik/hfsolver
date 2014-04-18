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

type umfpack_numeric
    type(c_ptr) :: numeric
end type

contains

subroutine factorize(n, Ap, Ai, Ax, umfpack_data)
integer, intent(in) :: n, Ap(:), Ai(:)
real(dp), intent(in) :: Ax(:)
type(umfpack_numeric), intent(out) :: umfpack_data
call stop_error("Not configured with UMFPACK.")
print *, n, Ap, Ai, Ax
umfpack_data%numeric = c_null_ptr
end subroutine

subroutine solve(Ap, Ai, Ax, x, b, umfpack_data)
integer, intent(in) :: Ap(:), Ai(:)
real(dp), intent(in) :: Ax(:), b(:)
real(dp), intent(out) :: x(:)
type(umfpack_numeric), intent(in) :: umfpack_data
type(c_ptr) :: tmp
call stop_error("Not configured with UMFPACK.")
print *, Ap, Ai, Ax, x, b
tmp = umfpack_data%numeric
end subroutine

subroutine free_data(umfpack_data)
type(umfpack_numeric), intent(inout) :: umfpack_data
call stop_error("Not configured with UMFPACK.")
umfpack_data%numeric = c_null_ptr
end subroutine

end module
