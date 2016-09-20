program test_ffte_par
use types, only: dp
use constants, only: pi
use fourier, only: dft, idft, fft, fft_vectorized, fft_pass, fft_pass_inplace, &
        fft_vectorized_inplace, calculate_factors, ifft_pass, fft2_inplace, &
        fft3_inplace, ifft3_inplace
use utils, only: assert, init_random, stop_error, get_int_arg, get_float_arg
use pffte, only: pfft3_init, pfft3, pifft3
use openmp, only: omp_get_wtime
use mpi2, only: mpi_finalize, MPI_COMM_WORLD, mpi_comm_rank, &
    mpi_comm_size, mpi_init, mpi_comm_split, MPI_INTEGER, &
    mpi_barrier, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_bcast, mpi_allreduce
use ofdft_fft, only: reciprocal_space_vectors, real_space_vectors, &
    real2fourier
implicit none

complex(dp), dimension(:,:,:), allocatable :: x3, x3d, x3d2, x3d3
real(dp), allocatable :: G2(:,:,:), X(:,:,:,:), X_global(:,:,:,:)
real(dp), allocatable :: G_global(:,:,:,:), G2_global(:,:,:)
real(dp) :: L(3), r2, alpha, Z, Eee, Eee_conv
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
        Ng = 32
        L = 2
    else
        if (command_argument_count() /= 9) then
            print *, "Usage:"
            print *,
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
!    print "(7i0)", myid, i, j, k, ijk_global
    X(i,j,k,:) = X_global(ijk_global(1), ijk_global(2), ijk_global(3), :)
    G2(i,j,k) = G2_global(ijk_global(1), ijk_global(2), ijk_global(3))
end do
end do
end do

! Setup two Gaussians with opposite charges, thus overall the density is net
! neutral:
alpha = 5
Z = 3
do k = 1, size(x3, 3)
do j = 1, size(x3, 2)
do i = 1, size(x3, 1)
    r2 =sum((X(i,j,k,:)-L/2)**2)
    x3(i, j, k) = Z*alpha**3/pi**(3._dp/2)*exp(-alpha**2*r2) &
        - Z*(alpha+1)**3/pi**(3._dp/2)*exp(-(alpha+1)**2*r2)
end do
end do
end do

Eee = pintegral(comm_all, L, real(x3, dp), Ng)
if (myid == 0) then
    print *, "integral:", myid, Eee
    call assert(abs(Eee) < 1e-10_dp)
end if
x3d = x3
call pfft3_init(Ng)
call mpi_barrier(comm_all, ierr)
call cpu_time(t1)
call pfft3(x3d, x3d2, commy, commz, Ng, nsub)
call mpi_barrier(comm_all, ierr)
!call real2fourier(x3d, x3d2)
call cpu_time(t2)
x3d2 = x3d2 / product(Ng)
!call pifft3(x3d2, x3d3, commy, commz, Ng, nsub)
!x3d3 = x3d3 * product(Ng)
call cpu_time(t3)
if (myid == 0) then
!    print *, myid, maxval(abs(x3 - x3d3))
    print *, "Timings (t, t*nproc):"
    print *, t2-t1, (t2-t1)*nproc
!    print *, t3-t2, (t3-t2)*nproc
end if
!call assert(all(abs(x3 - x3d3) < 5e-15_dp))

x3d3 = 4*pi*x3d2/G2
if (myid == 0) then
    x3d3(1,1,1) = 0
end if
Eee = pintegralG(comm_all, L, real(x3d3*conjg(x3d2), dp)) / 2 &
    - Z**2*alpha/sqrt(2*pi)
!Eee = pintegralG(comm_all, L, abs(x3d2)**2/G2)
if (myid == 0) then
    print *, "integralG:", myid, Eee
    Eee_conv = -17.465136801093962_dp
    print *, "error:", abs(Eee-Eee_conv)
    call assert(abs(Eee - Eee_conv) < 1e-12_dp)
end if
deallocate(x3, x3d, x3d2, x3d3)

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
