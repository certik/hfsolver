module mpi2
use types, only: dp
use mpi, only: MPI_COMM_WORLD, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_SUM, &
    mpi_comm_rank, mpi_comm_size, mpi_init, mpi_comm_split, mpi_barrier, &
    MPI_MAX, MPI_STATUS_SIZE
use mpi_dispatch, only: mpi_bcast_floats, mpi_bcast_float, mpi_bcast_ints, &
    mpi_bcast_int, mpi_allreduce_float, mpi_finalize => mpi_finalize_dispatch, &
    mpi_send_floats, mpi_recv_floats, mpi_send_float, mpi_recv_float
implicit none
private

! Implemented in the `mpi` module
public MPI_COMM_WORLD, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_MAX, &
    mpi_comm_rank, mpi_comm_size, mpi_init, mpi_comm_split, mpi_barrier, &
    MPI_STATUS_SIZE

! Implemented in the `mpi_dispatch` and interface blocks below
public mpi_bcast, mpi_finalize, mpi_allreduce, mpi_send, mpi_recv


interface mpi_bcast
    module procedure mpi_bcast_floats
    module procedure mpi_bcast_float
    module procedure mpi_bcast_ints
    module procedure mpi_bcast_int
end interface

interface mpi_allreduce
    module procedure mpi_allreduce_float
end interface

interface mpi_send
    module procedure mpi_send_float
    module procedure mpi_send_floats
end interface

interface mpi_recv
    module procedure mpi_recv_float
    module procedure mpi_recv_floats
end interface

end module
