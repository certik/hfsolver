module mpi_interface
use types, only: dp
use mpi
implicit none

interface
    subroutine mpi_finalize(ierr)
    integer, intent(out) :: ierr
    end subroutine
end interface

contains

subroutine mpi_bcast_floats(buffer, count_, datatype, root, comm, ierr)
integer, intent(in) :: count_
real(dp), intent(inout) :: buffer(count_)
integer, intent(in) :: datatype, root, comm
integer, intent(out) :: ierr
call mpi_bcast(buffer, count_, datatype, root, comm, ierr)
end subroutine

subroutine mpi_bcast_float(buffer, count_, datatype, root, comm, ierr)
integer, intent(in) :: count_
real(dp), intent(inout) :: buffer
integer, intent(in) :: datatype, root, comm
integer, intent(out) :: ierr
call mpi_bcast(buffer, count_, datatype, root, comm, ierr)
end subroutine

subroutine mpi_bcast_ints(buffer, count_, datatype, root, comm, ierr)
integer, intent(in) :: count_
integer, intent(inout) :: buffer(count_)
integer, intent(in) :: datatype, root, comm
integer, intent(out) :: ierr
call mpi_bcast(buffer, count_, datatype, root, comm, ierr)
end subroutine

subroutine mpi_bcast_int(buffer, count_, datatype, root, comm, ierr)
integer, intent(in) :: count_
integer, intent(inout) :: buffer
integer, intent(in) :: datatype, root, comm
integer, intent(out) :: ierr
call mpi_bcast(buffer, count_, datatype, root, comm, ierr)
end subroutine

end module
