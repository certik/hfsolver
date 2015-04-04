module mpi_interface
use types, only: dp
use mpi
implicit none

interface
    subroutine mpi_finalize(ierr)
    integer, intent(out) :: ierr
    end subroutine

    subroutine mpi_bcast(buffer, count_, datatype, root, comm, ierr)
    import :: dp
    integer, intent(in) :: count_
    real(dp), intent(inout) :: buffer(count_)
    integer, intent(in) :: datatype, root, comm
    integer, intent(out) :: ierr
    end subroutine
end interface

end module
