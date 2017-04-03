program test_parpack
use arpack, only: eig
use types, only: dp
use utils, only: stop_error, assert
use mpi2, only: mpi_finalize, MPI_COMM_WORLD, mpi_comm_rank, &
    mpi_comm_size, mpi_init, mpi_recv, mpi_send, MPI_DOUBLE_PRECISION, &
    MPI_STATUS_SIZE
use arpack, only: peig
implicit none

integer :: n, nloc, nev, ncv, i
real(dp), allocatable :: d(:), v(:,:)
real(dp), parameter :: e_ref(4) =  [5.4780178344413200_dp, &
    5.7004342714592289_dp, 5.8649444588087105_dp, 5.9659461993678038_dp]

integer :: comm, myid, nproc, ierr

call mpi_init(ierr)
comm  = MPI_COMM_WORLD
call mpi_comm_rank(comm, myid, ierr)
call mpi_comm_size(comm, nproc, ierr)

n = 16 ! global vector dimension
nev = 4
ncv = 8
nloc = n / nproc ! local vector dimension
call assert(nloc*nproc == n)
allocate(v(nloc,ncv), d(ncv))
call peig(comm, myid, nloc, nev, ncv, "LM", av, d, v)
do i = 1, nev
    if (myid == 0) print *, d(i), abs(d(i) - e_ref(i))
!    if (abs(d(i) - e_ref(i)) > 1e-14_dp) call stop_error("Error in ref. check")
end do

if (myid == 0) print *, "Done"

call mpi_finalize(ierr)

contains

  subroutine av(x, y)
  ! Compute y = A*x
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: y(:)
  real(dp) :: dd, dl, du, buf
  integer :: j, status(MPI_STATUS_SIZE)

  dd = 4
  dl = -1
  du = -1
  y(1) = dd*x(1) + du*x(2)
  do j = 2, size(y)-1
     y(j) = dl*x(j-1) + dd*x(j) + du*x(j+1)
  end do
  y(size(y)) = dl*x(size(y)-1) + dd*x(size(y))
  if (myid < nproc-1) then
     call mpi_send(x(size(x)), 1, MPI_DOUBLE_PRECISION, myid+1, myid+1, comm, &
         ierr)
  endif
  if (myid > 0 ) then
     call mpi_recv(buf, 1, MPI_DOUBLE_PRECISION, myid-1, myid, comm, &
         status, ierr)
     y(1) = y(1) + dl*buf
  endif
  if (myid > 0 ) then
     call mpi_send(x(1), 1, MPI_DOUBLE_PRECISION, myid-1, myid-1, comm, &
         ierr)
  endif
  if (myid < nproc-1) then
     call mpi_recv(buf, 1, MPI_DOUBLE_PRECISION, myid+1, myid, comm, &
         status, ierr)
     y(size(y)) = y(size(y)) + du*buf
  endif
  end

end
