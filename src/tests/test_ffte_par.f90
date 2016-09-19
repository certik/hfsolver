program test_ffte_par
use types, only: dp
use fourier, only: dft, idft, fft, fft_vectorized, fft_pass, fft_pass_inplace, &
        fft_vectorized_inplace, calculate_factors, ifft_pass, fft2_inplace, &
        fft3_inplace, ifft3_inplace
use utils, only: assert, init_random, stop_error, get_int_arg
use ffte, only: ffte_fft3_inplace => fft3_inplace, &
    ffte_ifft3_inplace => ifft3_inplace, pfft3, pfft3_init, pifft3
use openmp, only: omp_get_wtime
use mpi_interface, only: mpi_finalize, MPI_COMM_WORLD, mpi_comm_rank, &
    mpi_comm_size, mpi_init, mpi_comm_split, MPI_INTEGER, mpi_bcast_ints, &
    mpi_barrier
implicit none

complex(dp), dimension(:,:,:), allocatable :: x3, x3d, x3d2, x3d3
real(dp) :: tmp
integer :: l, m, n, i, j, k
integer :: Ng(3)
real(dp) :: t1, t2, t3

!  parallel variables
integer :: myid, comm_all, commy, commz, nproc, ierr, nsub(3), Ng_local(3)

call mpi_init(ierr)
comm_all  = MPI_COMM_WORLD
call mpi_comm_rank(comm_all, myid, ierr)
call mpi_comm_size(comm_all, nproc, ierr)
if (myid == 0) then
    if (command_argument_count() == 0) then
        nsub = 1
    else
        if (command_argument_count() /= 3) then
            print *, "Usage:"
            print *,
            print *, "test_ffte_par nsubx nsuby nsubz"
            call stop_error("Incorrect number of arguments.")
        end if
        nsub(1) = get_int_arg(1)
        nsub(2) = get_int_arg(2)
        nsub(3) = get_int_arg(3)
    end if
    print *, "nproc:", nproc
    print *, "nsub:", nsub

    if (product(nsub) /= nproc) then
        call stop_error("nproc must be equal to the number of subdomains")
    end if
    if (nsub(1) /= 1) then
        call stop_error("nsub(1) must be equal to 1")
    end if
end if
call mpi_bcast_ints(nsub, size(nsub), MPI_INTEGER, 0, comm_all, ierr)

call mpi_comm_split(comm_all, myid/nsub(2), 0, commy, ierr)
call mpi_comm_split(comm_all, mod(myid, nsub(2)), 0, commz, ierr)

Ng = 32
Ng_local = Ng / nsub

if (myid == 0) then
    print *, "Ng:      ", Ng
    print *, "Ng_local:", Ng_local
    if (.not. all(Ng_local * nsub == Ng)) then
        call stop_error("Ng must be equal to Ng_local * nsub")
    end if
end if
call mpi_barrier(comm_all, ierr)

allocate(x3(Ng_local(1), Ng_local(2), Ng_local(3)))
allocate(x3d(Ng_local(1), Ng_local(2), Ng_local(3)))
allocate(x3d2(Ng_local(1), Ng_local(2), Ng_local(3)))
allocate(x3d3(Ng_local(1), Ng_local(2), Ng_local(3)))
do k = 1, size(x3, 3)
do j = 1, size(x3, 2)
do i = 1, size(x3, 1)
    call random_number(tmp)
    x3(i, j, k) = tmp
end do
end do
end do
x3d = x3
call pfft3_init(Ng)
call cpu_time(t1)
call pfft3(x3d, x3d2, commy, commz, Ng, nsub)
call cpu_time(t2)
call pifft3(x3d2, x3d3, commy, commz, Ng, nsub)
call cpu_time(t3)
if (myid == 0) then
    print *, myid, maxval(abs(x3 - x3d3))
    print *, "Timings (t, t*nproc):"
    print *, t2-t1, (t2-t1)*nproc
    print *, t3-t2, (t3-t2)*nproc
end if
call assert(all(abs(x3 - x3d3) < 5e-15_dp))
deallocate(x3, x3d, x3d2, x3d3)

call mpi_finalize(ierr)

end program
