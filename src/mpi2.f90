module mpi2
use types, only: dp
use mpi
use mpi_dispatch, only: mpi_bcast_floats, mpi_bcast_float, mpi_bcast_ints, &
    mpi_bcast_int, mpi_allreduce_float
implicit none
private
public mpi_bcast, mpi_finalize, MPI_COMM_WORLD, mpi_comm_rank, mpi_comm_size, &
    mpi_init, mpi_comm_split, MPI_INTEGER, mpi_barrier, MPI_DOUBLE_PRECISION, &
    MPI_SUM, mpi_allreduce

interface
    subroutine mpi_finalize(ierr)
    integer, intent(out) :: ierr
    end subroutine
end interface

interface mpi_bcast
    module procedure mpi_bcast_floats
    module procedure mpi_bcast_float
    module procedure mpi_bcast_ints
    module procedure mpi_bcast_int
end interface

interface mpi_allreduce
    module procedure mpi_allreduce_float
end interface

end module
