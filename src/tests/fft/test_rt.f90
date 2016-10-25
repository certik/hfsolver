program test_rt
use types, only: dp
use constants, only: Ha2eV, i_
use fourier, only: dft, idft, fft, fft_vectorized, fft_pass, fft_pass_inplace, &
        fft_vectorized_inplace, calculate_factors, ifft_pass, fft2_inplace, &
        fft3_inplace, ifft3_inplace
use utils, only: assert, init_random, stop_error, get_int_arg, get_float_arg, &
    allocate_mold, clock
use ffte, only: factor
use ofdft, only: read_pseudo
use pofdft_fft, only: pfft3_init, preal2fourier, pfourier2real, &
    real_space_vectors, reciprocal_space_vectors, calculate_myxyz, &
    pintegral, pintegralG, free_energy, free_energy_min, &
    radial_potential_fourier, psum, pmaxval
use openmp, only: omp_get_wtime
use mpi2, only: mpi_finalize, MPI_COMM_WORLD, mpi_comm_rank, &
    mpi_comm_size, mpi_init, mpi_comm_split, MPI_INTEGER, &
    mpi_barrier, mpi_bcast
use md, only: positions_bcc
implicit none

complex(dp), dimension(:,:,:), allocatable :: neG, VenG, psiG, psi, tmp
real(dp), dimension(:,:,:), allocatable :: G2, Hn, Ven0G, ne, Ven
real(dp), allocatable :: G(:,:,:,:), X(:,:,:,:), Xion(:,:), R(:), Ven_rad(:), &
    current(:,:,:,:)
complex(dp), allocatable :: dpsi(:,:,:,:)
real(dp) :: L(3), Z
integer :: i, j
integer :: Ng(3)
integer :: LNPU(3)
integer :: cg_iter, natom, u
real(dp) :: T_eV, T_au, Eee, Een, Ts, Exc, Etot, Ediff, V0, mu, &
    mu_Hn, dt, psi_norm, t, Hn_mu_diff, EvW, &
    lambda, current_avg(3), A0, lambdaK

!  parallel variables
integer :: comm_all, commy, commz, nproc, ierr, nsub(3), Ng_local(3)
integer :: myid ! my ID (MPI rank), starts from 0
integer :: myxyz(3) ! myid, converted to the (x, y, z) box, starts from 0

T_eV = 34.5_dp
T_au = T_ev / Ha2eV
natom = 128
L = 8.1049178668765851_dp
lambda = 0
lambdaK = 1 ! Coefficient of the Laplace operator in the kinetic term

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
    else
        if (command_argument_count() /= 6) then
            print *, "Usage:"
            print *
            print *, "test_ffte_par Ng(3) nsub(3)"
            call stop_error("Incorrect number of arguments.")
        end if
        Ng(1) = get_int_arg(1)
        Ng(2) = get_int_arg(2)
        Ng(3) = get_int_arg(3)
        nsub(1) = get_int_arg(4)
        nsub(2) = get_int_arg(5)
        nsub(3) = get_int_arg(6)
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
end if
call mpi_bcast(nsub, size(nsub), MPI_INTEGER, 0, comm_all, ierr)
call mpi_bcast(Ng, size(Ng), MPI_INTEGER, 0, comm_all, ierr)
call mpi_bcast(Ng_local, size(Ng_local), MPI_INTEGER, 0, comm_all, ierr)
call pfft3_init(myid, comm_all, Ng, nsub)

myxyz = calculate_myxyz(myid, nsub)

! Note that myxyz(3) corresponds to commy, and myxyz(2) to commz
call mpi_comm_split(comm_all, myxyz(3), 0, commy, ierr)
call mpi_comm_split(comm_all, myxyz(2), 0, commz, ierr)


allocate(ne(Ng_local(1), Ng_local(2), Ng_local(3)))
allocate(neG(Ng_local(1), Ng_local(2), Ng_local(3)))
call allocate_mold(G2, ne)
call allocate_mold(Hn, ne)
call allocate_mold(Ven, ne)
call allocate_mold(Ven0G, ne)
call allocate_mold(VenG, neG)
call allocate_mold(psiG, neG)
call allocate_mold(psi, neG)
call allocate_mold(tmp, neG)
allocate(X(Ng_local(1), Ng_local(2), Ng_local(3), 3))
allocate(dpsi(Ng_local(1), Ng_local(2), Ng_local(3), 3))
call allocate_mold(G, X)
call allocate_mold(current, X)
if (myid == 0) print *, "Load initial position"
allocate(Xion(3, natom))
! For now assume a box, until positions_bcc can accept a vector L(:)
! And radial_potential_fourier
call assert(abs(L(2)-L(1)) < 1e-15_dp)
call assert(abs(L(3)-L(1)) < 1e-15_dp)
call positions_bcc(Xion, L(1))
call real_space_vectors(L, X, Ng, myxyz)
call reciprocal_space_vectors(L, G, G2, Ng, myxyz)
if (myid == 0) print *, "Radial nuclear potential FFT"
call read_pseudo("../fem/D.pseudo", R, Ven_rad, Z, Ediff)
call radial_potential_fourier(R, Ven_rad, L, Z, G2, Ven0G, V0)
if (myid == 0) print *, "    Done."

VenG = 0
do i = 1, natom
    VenG = VenG - Ven0G * exp(-i_ * &
        (G(:,:,:,1)*Xion(1,i) + G(:,:,:,2)*Xion(2,i) + G(:,:,:,3)*Xion(3,i)))
end do


! Do CG minimization
if (myid == 0) print *, "CG minimization:"
ne = natom / product(L)
call free_energy_min(myid, comm_all, commy, commz, Ng, nsub, &
        real(natom, dp), natom, L, G2, T_au, VenG, ne, 1e-12_dp, &
        Eee, Een, Ts, Exc, Etot, cg_iter, .true., lambda, EvW)
call free_energy(myid, comm_all, commy, commz, Ng, nsub, &
        L, G2, T_au, VenG, ne, Eee, Een, Ts, Exc, Etot, Hn, &
        .true., .true., .true., lambda, lambda, EvW)

mu = 1._dp / natom * pintegral(comm_all, L, ne * Hn, Ng)
mu_Hn = psum(comm_all, Hn)/product(Ng)
if (myid == 0) then
    print *, mu, mu_Hn
    print *, "Summary of energies [a.u.]:"
    print "('    EvW  = ', f14.8)", EvW
    print "('    Ts   = ', f14.8)", Ts
    print "('    Een  = ', f14.8)", Een
    print "('    Eee  = ', f14.8)", Eee
    print "('    Exc  = ', f14.8)", Exc
    print *, "   ---------------------"
    print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV
end if


call free_energy(myid, comm_all, commy, commz, Ng, nsub, &
        L, G2, T_au, VenG, ne, Eee, Een, Ts, Exc, Etot, Hn, &
        .true., .true., .true., lambda, lambda-lambdaK, EvW)

mu = psum(comm_all, Hn)/product(Ng)
Hn_mu_diff = pmaxval(comm_all, abs(Hn - mu))
if (myid == 0) then
    print *, "Etot =", Etot
    print *, "mu = ", mu
    print *, "max(abs(H-mu)) = ", Hn_mu_diff
end if
call assert(all(ne > 0))

if (myid == 0) then
    print *
    print *, "------------------------------------------------------------------"
    print *, "Propagation"
end if

! Propagate

if (myid == 0) then
    print *, "E_max =", maxval(abs(Hn)), "; dt <", 1/maxval(abs(Hn))
end if

dt = 1e-3_dp

if (myid == 0) then
    print *, "dt =", dt
end if

! Do first step by hand:
if (myid == 0) print *, "First step"
psi = sqrt(ne)

t = 0
A0 = 1e-3_dp

psi = psi * exp(-i_*Hn*dt/2)
call preal2fourier(psi, psiG, commy, commz, Ng, nsub)
psiG = psiG * exp(-i_*G2*dt/2*lambdaK)
call pfourier2real(psiG, psi, commy, commz, Ng, nsub)
psi = psi * exp(-i_*Hn*dt/2)

ne = abs(psi)**2
psi_norm = pintegral(comm_all, L, ne, Ng)
if (myid == 0) print *, "norm of psi:", psi_norm

if (myid == 0) open(newunit=u, file="of_cond.txt", status="replace")
do i = 1, 10
    t = t + dt
    if (myid == 0) print *, "iter =", i, "time =", t
    psi = psi * exp(-i_*Hn*dt/2)
    call preal2fourier(psi, psiG, commy, commz, Ng, nsub)
    psiG = psiG * exp(-i_*G2*dt/2*lambdaK)
    psiG = psiG * exp(A0*G(:,:,:,1)*dt)
    call pfourier2real(psiG, psi, commy, commz, Ng, nsub)
    psi = psi * exp(-i_*Hn*dt/2)
    ne = abs(psi)**2
    psi_norm = pintegral(comm_all, L, ne, Ng)
    if (myid == 0) print *, "norm of psi:", psi_norm

    call free_energy(myid, comm_all, commy, commz, Ng, nsub, &
            L, G2, T_au, VenG, ne, Eee, Een, Ts, Exc, Etot, Hn, &
            .true., .true., .true., lambda, lambda-lambdaK, EvW)

    !Calculate Current
    call preal2fourier(psi, psiG, commy, commz, Ng, nsub)
    do j = 1, 3
        call pfourier2real(i_*G(:,:,:,j)*psiG, dpsi(:,:,:,j), &
                commy, commz, Ng, nsub)
        tmp = (conjg(psi)*dpsi(:,:,:,j)-psi*conjg(dpsi(:,:,:,j))) / (2*natom*i_)
        if (maxval(abs(aimag(tmp))) > 1e-12_dp) then
            print *, "INFO: current  max imaginary part:", maxval(aimag(tmp))
        end if
        current(:,:,:,j) = real(tmp, dp)
        current_avg(j) = pintegral(comm_all, L, current(:,:,:,j), Ng)
    end do


    mu = 1._dp / natom * pintegral(comm_all, L, ne * Hn, Ng)
    if (myid == 0) then
        print *, mu
        print *, "Summary of energies [a.u.]:"
        print "('    EvW  = ', f14.8)", EvW
        print "('    Ts   = ', f14.8)", Ts
        print "('    Een  = ', f14.8)", Een
        print "('    Eee  = ', f14.8)", Eee
        print "('    Exc  = ', f14.8)", Exc
        print *, "   ---------------------"
        print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV
        write(u,*) i, t, current_avg
    end if
end do
if (myid == 0) close(u)
if (myid == 0) print *, "Done"

call mpi_finalize(ierr)
end program
