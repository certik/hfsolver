module mpi_dispatch
! This module contains subroutines that are not defined in the `mpi` module.
! They are used implicitly here, because it is not possible to create an
! interface declaration that correctly dispatches on the type. So this module
! will give "implicit interface" compiler warning. But the rest of the code
! uses properly typed and declared subroutines, and they have the same names as
! in the MPI 2 standard, thanks to the mpi2.f90 interface blocks. Users should
! just import the mpi2 module.
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

! mpi_finalize is available in openmpi, but not mpich
subroutine mpi_finalize_dispatch(ierr)
integer, intent(out) :: ierr
call mpi_finalize(ierr)
end subroutine

! mpi_send

subroutine mpi_send_float(buffer, n, datatype, dest, tag, comm, ierr)
integer, intent(in) :: n
real(dp), intent(in) :: buffer
integer, intent(in) :: datatype, dest, tag, comm
integer, intent(out) :: ierr
call mpi_send(buffer, n, datatype, dest, tag, comm, ierr)
end subroutine

subroutine mpi_send_floats(buffer, n1, n2, n3, datatype, dest, tag, comm, ierr)
integer, intent(in) :: n1, n2, n3
real(dp), intent(in) :: buffer(n1,n2,n3)
integer, intent(in) :: datatype, dest, tag, comm
integer, intent(out) :: ierr
call mpi_send(buffer, n1*n2*n3, datatype, dest, tag, comm, ierr)
end subroutine

! mpi_recv

subroutine mpi_recv_float(buffer, n, datatype, source, tag, comm, &
        stat, ierr)
use mpi, only: MPI_STATUS_SIZE
integer, intent(in) :: n
real(dp), intent(out) :: buffer
integer, intent(in) :: datatype, source, tag, comm
integer, intent(out) :: stat(MPI_STATUS_SIZE), ierr
call mpi_recv(buffer, n, datatype, source, tag, comm, stat, ierr)
end subroutine

subroutine mpi_recv_floats(buffer, n1, n2, n3, datatype, source, tag, comm, &
        stat, ierr)
use mpi, only: MPI_STATUS_SIZE
integer, intent(in) :: n1, n2, n3
real(dp), intent(out) :: buffer(n1,n2,n3)
integer, intent(in) :: datatype, source, tag, comm
integer, intent(out) :: stat(MPI_STATUS_SIZE), ierr
call mpi_recv(buffer, n1*n2*n3, datatype, source, tag, comm, stat, ierr)
end subroutine

end module
