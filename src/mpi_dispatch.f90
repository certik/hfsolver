module mpi_dispatch
use types, only: dp
implicit none

contains

! mpi_bcast

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

! mpi_allreduce

subroutine mpi_allreduce_float(sendbuf, recvbuf, count_, datatype, op, comm, ierr)
integer, intent(in) :: count_
real(dp), intent(in) :: sendbuf
real(dp), intent(out) :: recvbuf
integer, intent(in) :: datatype, op, comm
integer, intent(out) :: ierr
call mpi_allreduce(sendbuf, recvbuf, count_, datatype, op, comm, ierr)
end subroutine

end module
