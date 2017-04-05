program test_psch3d
use types, only: dp
use constants, only: Ha2eV, i_, density2gcm3, u2au, s2au, K2au
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
    radial_potential_fourier, psum, pmaxval, collate
use openmp, only: omp_get_wtime
use mpi2, only: mpi_finalize, MPI_COMM_WORLD, mpi_comm_rank, &
    mpi_comm_size, mpi_init, mpi_comm_split, MPI_INTEGER, &
    mpi_barrier, mpi_bcast
use md, only: positions_bcc, positions_fcc
use arpack, only: peig, eig
use pksdft_fft, only: solve_schroedinger
implicit none

complex(dp), dimension(:,:,:), allocatable :: neG, VenG, psiG, psi, tmp
real(dp), dimension(:,:,:), allocatable :: G2, Hn, Htot, HtotG, Ven0G, ne, Ven
real(dp), allocatable :: G(:,:,:,:), X(:,:,:,:), Xion(:,:), R(:), Ven_rad(:), &
    current(:,:,:,:), tmp_global(:,:,:), eigs(:), orbitals(:,:,:,:)
complex(dp), allocatable :: dpsi(:,:,:,:)
real(dp) :: L(3), Z, omega
integer :: i
integer :: Ng(3)
integer :: LNPU(3)
integer :: natom, u2
integer :: Lmul
logical :: velocity_gauge
real(dp) :: T_eV, T_au, Etot, Ediff, V0, mu, &
    dt, psi_norm, t, &
    alpha, rho, Ekin, Epot, Etot_last
real(dp), allocatable :: m(:)
real(dp) :: t1, t2, eps
integer :: nev, ncv

!  parallel variables
integer :: comm_all, commy, commz, nproc, ierr, nsub(3), Ng_local(3)
integer :: myid ! my ID (MPI rank), starts from 0
integer :: myxyz(3) ! myid, converted to the (x, y, z) box, starts from 0

rho = 0.01_dp / density2gcm3  ! g/cc
T_eV = 50._dp
T_au = T_ev / Ha2eV
natom = 1
dt = 1e-4_dp
alpha = 137
allocate(m(natom))
m = 2._dp * u2au ! Using Argon mass in atomic mass units [u]
velocity_gauge = .true. ! velocity or length gauge?

L = (sum(m) / rho)**(1._dp/3)
Lmul = 1
L = L * Lmul
rho = sum(m) / product(L)

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
        Ng = 4 * Lmul
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

    print *, "Input:"
    print *, "N =", natom
    print *, "rho = ", rho * density2gcm3, "g/cm^3 = ", rho, "a.u."
    print *, "T =", T_au / K2au, "K =", T_au, "a.u. =", T_au * Ha2eV, "eV"
    print "('dt =', es10.2, ' s = ', es10.2, ' ps = ', es10.2, ' a.u.')", &
        dt/s2au, dt/s2au * 1e12_dp, dt
    print *
    print *, "Calculated quantities:"
    print *, "L =", L, "a.u."
    print *
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


allocate(tmp_global(Ng(1), Ng(2), Ng(3)))
allocate(ne(Ng_local(1), Ng_local(2), Ng_local(3)))
allocate(neG(Ng_local(1), Ng_local(2), Ng_local(3)))
call allocate_mold(G2, ne)
call allocate_mold(Hn, ne)
call allocate_mold(Htot, ne)
call allocate_mold(HtotG, ne)
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
!call positions_fcc(Xion, L(1))
Xion(:, 1) = L/2
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

call pfourier2real(VenG, Ven, commy, commz, Ng, nsub)
call collate(comm_all, myid, nsub, 0, Ven, tmp_global)
if (myid == 0) then
    open(newunit=u2, file="sch_pot.txt", status="replace")
    write(u2,*) tmp_global
    close(u2)
end if
omega = 1.123_dp
!do k = 1, Ng_local(3)
!do j = 1, Ng_local(2)
!do i = 1, Ng_local(1)
!    Ven(i,j,k) = omega**2*(sqrt(sum((X(i,j,k,:)-L/2)**2)))**2 / 2
!end do
!end do
!end do
Hn = Ven

if (myid == 0) print *, "IT minimization:"
ne = natom / product(L)
psi = sqrt(ne)

t = 0
dt = 2e-1_dp

eps = 1e-20_dp
Etot_last = -1e9_dp
call cpu_time(t1)
i = 0
do
    i = i + 1
    t = t + dt
!    if (myid == 0) print *, "iter =", i, "time =", t
    psi = psi * exp(-(Hn)*dt/2)
    call preal2fourier(psi, psiG, commy, commz, Ng, nsub)
    psiG = psiG * exp(-G2*dt/2)
    call pfourier2real(psiG, psi, commy, commz, Ng, nsub)
    psi = psi * exp(-(Hn)*dt/2)
    ne = abs(psi)**2
    psi_norm = pintegral(comm_all, L, ne, Ng)
!    if (myid == 0) print *, "Initial norm of psi:", psi_norm
    psi = sqrt(natom / psi_norm) * psi
    ne = abs(psi)**2
    psi_norm = pintegral(comm_all, L, ne, Ng)
!    if (myid == 0) print *, "norm of psi:", psi_norm

    mu = 1._dp / natom * pintegral(comm_all, L, ne * Hn, Ng)
    call preal2fourier(psi, psiG, commy, commz, Ng, nsub)
    Ekin = 1._dp/2 * pintegralG(comm_all, L, G2*abs(psiG)**2)
    Epot = pintegral(comm_all, L, Hn*ne, Ng)
    Etot = Ekin + Epot
    if (myid == 0) then
!        print *, mu
!        print *, "Summary of energies [a.u.]:"
!        print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV
!        print *, "Exact:", 3*omega/2
        print *, i, Etot, dt
    end if
    if (abs(Etot_last - Etot) < eps .and. i > 3) then
        dt = dt/10
        i = 1
    end if
    Etot_last = Etot
    if (i > 20000 .or. dt < 1e-12_dp) then
        exit
    end if
end do
call cpu_time(t2)
if (myid == 0) print *, "IT Time:", t2-t1

if (myid == 0) print *, "Arpack"

nev = 4
ncv = 64
allocate(eigs(nev), orbitals(Ng_local(1),Ng_local(2),Ng_local(3),nev))
call solve_schroedinger(myid, comm_all, commy, commz, Ng, nsub, Ven, &
        L, G2, nev, ncv, eigs, orbitals)
if (myid == 0) then
    print *, "Eigenvalues:"
    do i = 1, nev
        print *, i, eigs(i)
    end do
    print *, "Arpack vs IT:"
    print *, abs(eigs(1) - Etot)
    call assert(abs(eigs(1) - Etot) < 1e-12_dp)
end if

if (myid == 0) print *, "Done"

call mpi_finalize(ierr)

end program
