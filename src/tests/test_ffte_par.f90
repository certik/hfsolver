program test_ffte_par
use types, only: dp
use constants, only: pi
use fourier, only: dft, idft, fft, fft_vectorized, fft_pass, fft_pass_inplace, &
        fft_vectorized_inplace, calculate_factors, ifft_pass, fft2_inplace, &
        fft3_inplace, ifft3_inplace
use utils, only: assert, init_random, stop_error, get_int_arg, get_float_arg
use ffte, only: factor
use pofdft_fft, only: pfft3_init, preal2fourier, pfourier2real, &
    real_space_vectors, reciprocal_space_vectors
use openmp, only: omp_get_wtime
use mpi2, only: mpi_finalize, MPI_COMM_WORLD, mpi_comm_rank, &
    mpi_comm_size, mpi_init, mpi_comm_split, MPI_INTEGER, &
    mpi_barrier, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_bcast, mpi_allreduce
implicit none

complex(dp), dimension(:,:,:), allocatable :: ne, neG, ne2
real(dp), allocatable :: G(:,:,:,:), G2(:,:,:), X(:,:,:,:)
real(dp) :: L(3), r2, alpha, Z, Eee, Eee_conv
integer :: i, j, k
integer :: Ng(3)
real(dp) :: t1, t2, t3
integer :: LNPU(3)

!  parallel variables
integer :: comm_all, commy, commz, nproc, ierr, nsub(3), Ng_local(3)
integer :: myid ! my ID (MPI rank), starts from 0
integer :: myxyz(3) ! myid, converted to the (x, y, z) box, starts from 0

call mpi_init(ierr)
comm_all  = MPI_COMM_WORLD
call mpi_comm_rank(comm_all, myid, ierr)
call mpi_comm_size(comm_all, nproc, ierr)
if (myid == 0) then
    if (command_argument_count() == 0) then
        call factor(nproc, LNPU)
        nsub(3) = (2**(LNPU(1)/2))*(3**(LNPU(2)/2))*(5**(LNPU(3)/2))
        nsub(2) = nproc / nsub(3)
        nsub(1) = 1
        Ng = 32
        L = 2
    else
        if (command_argument_count() /= 9) then
            print *, "Usage:"
            print *
            print *, "test_ffte_par L(3) Ng(3) nsub(3)"
            call stop_error("Incorrect number of arguments.")
        end if
        L(1) = get_float_arg(1)
        L(2) = get_float_arg(2)
        L(3) = get_float_arg(3)
        Ng(1) = get_int_arg(4)
        Ng(2) = get_int_arg(5)
        Ng(3) = get_int_arg(6)
        nsub(1) = get_int_arg(7)
        nsub(2) = get_int_arg(8)
        nsub(3) = get_int_arg(9)
    end if
    Ng_local = Ng / nsub

    print *, "L:       ", L
    print *, "nproc:   ", nproc
    print *, "nsub:    ", nsub
    print *, "Ng:      ", Ng
    print *, "Ng_local:", Ng_local

    if (product(nsub) /= nproc) then
        call stop_error("nproc must be equal to the number of subdomains")
    end if
    if (nsub(1) /= 1) then
        call stop_error("nsub(1) must be equal to 1")
    end if
    if (.not. all(Ng_local * nsub == Ng)) then
        call stop_error("Ng must be equal to Ng_local * nsub")
    end if
end if
call mpi_bcast(L, size(L), MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
call mpi_bcast(nsub, size(nsub), MPI_INTEGER, 0, comm_all, ierr)
call mpi_bcast(Ng, size(Ng), MPI_INTEGER, 0, comm_all, ierr)
call mpi_bcast(Ng_local, size(Ng_local), MPI_INTEGER, 0, comm_all, ierr)

myxyz = [0, mod(myid, nsub(2)), myid/nsub(2)]

call mpi_comm_split(comm_all, myxyz(3), 0, commy, ierr)
call mpi_comm_split(comm_all, myxyz(2), 0, commz, ierr)


allocate(ne(Ng_local(1), Ng_local(2), Ng_local(3)))
allocate(neG(Ng_local(1), Ng_local(2), Ng_local(3)))
allocate(ne2(Ng_local(1), Ng_local(2), Ng_local(3)))
allocate(X(Ng_local(1), Ng_local(2), Ng_local(3), 3))
allocate(G(Ng_local(1), Ng_local(2), Ng_local(3), 3))
allocate(G2(Ng_local(1), Ng_local(2), Ng_local(3)))
call real_space_vectors(L, X, Ng, myxyz)
call reciprocal_space_vectors(L, G, G2, Ng, myxyz)
G(1, 1, 1, :) = 1 ! To avoid division by 0
G2(1, 1, 1) = 1 ! To avoid division by 0

! Setup two Gaussians with opposite charges, thus overall the density is net
! neutral:
alpha = 5
Z = 3
do k = 1, size(ne, 3)
do j = 1, size(ne, 2)
do i = 1, size(ne, 1)
    r2 =sum((X(i,j,k,:)-L/2)**2)
    ne(i, j, k) = Z*alpha**3/pi**(3._dp/2)*exp(-alpha**2*r2) &
        - Z*(alpha+1)**3/pi**(3._dp/2)*exp(-(alpha+1)**2*r2)
end do
end do
end do

Eee = pintegral(comm_all, L, real(ne, dp), Ng)
if (myid == 0) then
    print *, "integral:", myid, Eee
    call assert(abs(Eee) < 1e-10_dp)
end if
ne2 = 0
call pfft3_init(Ng)
call mpi_barrier(comm_all, ierr)
call cpu_time(t1)
call preal2fourier(ne, neG, commy, commz, Ng, nsub)
call mpi_barrier(comm_all, ierr)
call cpu_time(t2)
call pfourier2real(neG, ne2, commy, commz, Ng, nsub)
call mpi_barrier(comm_all, ierr)
call cpu_time(t3)
if (myid == 0) then
    print *, "Timings (t, t*nproc):"
    print *, t2-t1, (t2-t1)*nproc
    print *, t3-t2, (t3-t2)*nproc
    print *, "Error:", maxval(abs(ne-ne2))
    call assert(all(abs(ne - ne2) < 1e-12_dp))
end if
call mpi_barrier(comm_all, ierr)

if (myid == 0) then
    neG(1,1,1) = 0
end if
Eee = pintegralG(comm_all, L, 2*pi*abs(neG)**2/G2) - Z**2*alpha/sqrt(2*pi)
if (myid == 0) then
    print *, "integralG:", myid, Eee
    Eee_conv = -17.465136801093962_dp
    print *, "error:", abs(Eee-Eee_conv)
    call assert(abs(Eee - Eee_conv) < 1e-12_dp)
end if

call mpi_finalize(ierr)

contains

real(dp) function pintegral(comm, L, f, Ng) result(r)
! Calculates the integral over 'f' in parallel, returns the answer on all
! processors.
integer, intent(in) :: comm
real(dp), intent(in) :: L(:), f(:,:,:)
integer, intent(in) :: Ng(:)
real(dp) :: myr
myr = sum(f)
call mpi_allreduce(myr, r, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
r = r * product(L/Ng)
end function

real(dp) function pintegralG(comm, L, fG) result(r)
! Calculates the integral over 'fG' in reciprocal space in parallel, returns
! the answer on all processors.
integer, intent(in) :: comm
real(dp), intent(in) :: L(:)
real(dp), intent(in) :: fG(:, :, :)
real(dp) :: myr
!myr = sum(real(fG, dp))
myr = sum(fG)
call mpi_allreduce(myr, r, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr)
r = r * product(L)
end function

end program
