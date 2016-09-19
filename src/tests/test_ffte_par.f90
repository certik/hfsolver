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
    mpi_barrier, MPI_DOUBLE_PRECISION, MPI_SUM
use ofdft_fft, only: reciprocal_space_vectors, real_space_vectors
implicit none

complex(dp), dimension(:,:,:), allocatable :: x3, x3d, x3d2, x3d3
real(dp), allocatable :: G2(:,:,:), X(:,:,:,:), X_global(:,:,:,:)
real(dp), allocatable :: G_global(:,:,:,:), G2_global(:,:,:)
real(dp) :: L, r2
integer :: i, j, k
integer :: Ng(3)
real(dp) :: t1, t2, t3
integer :: LNPU(3)

!  parallel variables
integer :: comm_all, commy, commz, nproc, ierr, nsub(3), Ng_local(3)
integer :: myid ! my ID (MPI rank), starts from 0
integer :: myxyz(3) ! myid, converted to the (x, y, z) box, starts from 0
integer :: ijk_global(3)

call mpi_init(ierr)
comm_all  = MPI_COMM_WORLD
call mpi_comm_rank(comm_all, myid, ierr)
call mpi_comm_size(comm_all, nproc, ierr)
if (myid == 0) then
    if (command_argument_count() == 0) then
        call FACTOR(nproc, LNPU)
        nsub(3) = (2**(LNPU(1)/2))*(3**(LNPU(2)/2))*(5**(LNPU(3)/2))
        nsub(2) = nproc / nsub(3)
        nsub(1) = 1
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

myxyz = [0, myid/nsub(2), mod(myid, nsub(2))]

call mpi_comm_split(comm_all, myxyz(2), 0, commy, ierr)
call mpi_comm_split(comm_all, myxyz(3), 0, commz, ierr)

myxyz = [0, mod(myid, nsub(2)), myid/nsub(2)]

Ng = 32
Ng_local = Ng / nsub
L = 10

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
allocate(G_global(Ng(1), Ng(2), Ng(3), 3))
allocate(G2_global(Ng(1), Ng(2), Ng(3)))
allocate(X_global(Ng(1), Ng(2), Ng(3), 3))
allocate(G2(Ng_local(1), Ng_local(2), Ng_local(3)))
allocate(X(Ng_local(1), Ng_local(2), Ng_local(3), 3))
call real_space_vectors(L, X_global)
call reciprocal_space_vectors(L, G_global, G2_global)
! Convert X_global to X (local)
! Convert G2_global to G2 (local)
do k = 1, size(x3, 3)
do j = 1, size(x3, 2)
do i = 1, size(x3, 1)
    ijk_global = [i, j, k] + myxyz*Ng_local
    X(i,j,k,:) = X_global(ijk_global(1), ijk_global(2), ijk_global(3), :)
    G2(i,j,k) = G2_global(ijk_global(1), ijk_global(2), ijk_global(3))
end do
end do
end do

do k = 1, size(x3, 3)
do j = 1, size(x3, 2)
do i = 1, size(x3, 1)
    r2 =sum((X(i,j,k,:)-[L/2,L/2,L/2])**2)
    x3(i, j, k) = exp(-r2)
end do
end do
end do

print *, "integral:", myid, pintegral(comm_all, L, real(x3, dp), Ng)
x3d = x3
call pfft3_init(Ng)
call cpu_time(t1)
call pfft3(x3d, x3d2, commy, commz, Ng, nsub)
x3d2 = x3d2 / product(Ng)
call cpu_time(t2)
call pifft3(x3d2, x3d3, commy, commz, Ng, nsub)
x3d3 = x3d3 * product(Ng)
call cpu_time(t3)
if (myid == 0) then
    print *, myid, maxval(abs(x3 - x3d3))
    print *, "Timings (t, t*nproc):"
    print *, t2-t1, (t2-t1)*nproc
    print *, t3-t2, (t3-t2)*nproc
end if
call assert(all(abs(x3 - x3d3) < 5e-15_dp))

if (myid == 0) then
    x3d2(1,1,1) = 0
end if
x3d2 = x3d2/G2 ! FIXME: this depends on the decomposition...
! Probably x3d2 is distributed differently than G2 for some reason.
print *, "integralG:", myid, pintegralG(comm_all, L, x3d2)
deallocate(x3, x3d, x3d2, x3d3)

call mpi_finalize(ierr)

contains

real(dp) function pintegral(comm, L, f, Ng) result(r)
! Calculates the integral over 'f' in parallel, returns the answer on all
! processors.
integer, intent(in) :: comm
real(dp), intent(in) :: L, f(:,:,:)
integer, intent(in) :: Ng(:)
real(dp) :: myr
myr = sum(f)
call mpi_allreduce(myr, r, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
r = r * L**3 / product(Ng)
end function

real(dp) function pintegralG(comm, L, fG) result(r)
! Calculates the integral over 'fG' in reciprocal space in parallel, returns
! the answer on all processors.
integer, intent(in) :: comm
real(dp), intent(in) :: L
complex(dp), intent(in) :: fG(:, :, :)
real(dp) :: myr
myr = sum(real(fG, dp))
call mpi_allreduce(myr, r, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
r = r * L**3
end function

end program
