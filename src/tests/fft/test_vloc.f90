program test_vloc
use types, only: dp
use constants, only: Ha2eV, density2gcm3, u2au, s2au, K2au, i_, pi
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
    radial_potential_fourier, psum, pmaxval, collate, poisson_kernel
use openmp, only: omp_get_wtime
use mpi2, only: mpi_finalize, MPI_COMM_WORLD, mpi_comm_rank, &
    mpi_comm_size, mpi_init, mpi_comm_split, MPI_INTEGER, &
    mpi_barrier, mpi_bcast
use md, only: positions_bcc, positions_fcc
use arpack, only: peig, eig
use pksdft_fft, only: solve_schroedinger
use xc, only: xc_pz
use mixings, only: mixing_linear, mixing_linear_adapt
use ewald_sums, only: ewald_box
implicit none

complex(dp), dimension(:,:,:), allocatable :: neG, psiG, psi, tmp
real(dp), dimension(:,:,:), allocatable :: G2, Htot, HtotG, Ven0G, ne, Vloc, &
    Veff
real(dp), allocatable :: G(:,:,:,:), X(:,:,:,:), Xion(:,:), q(:), &
    current(:,:,:,:), eigs(:), orbitals(:,:,:,:), eigs_ref(:), occ(:), &
    Vee(:,:,:), Vee_xc(:,:,:), exc(:,:,:), Vxc(:,:,:), forces(:,:), &
    cutfn(:,:,:)
complex(dp), allocatable :: dpsi(:,:,:,:), VeeG(:,:,:), VenG(:,:,:)
real(dp) :: L(3), r, stress(6)
integer :: i, j, k, u
integer :: Ng(3)
integer :: LNPU(3)
integer :: natom
logical :: velocity_gauge
real(dp) :: T_eV, T_au, dt, alpha, rho, norm, w2, Vmin, Ekin, Etot, &
    Eee, Een_loc, E_xc, Enn, Een_core
real(dp) :: rloc, C1, C2, Zion
real(dp), allocatable :: m(:)
integer :: nev, ncv, na
real(dp), parameter :: D(5) = [0.65435_dp, 2.45106_dp, -1.536643785333E-01_dp, &
    1.153664378533E+00_dp, 5.0000_dp]

!  parallel variables
integer :: comm_all, commy, commz, nproc, ierr, nsub(3), Ng_local(3)
integer :: myid ! my ID (MPI rank), starts from 0
integer :: myxyz(3) ! myid, converted to the (x, y, z) box, starts from 0

rho = 0.01_dp / density2gcm3  ! g/cc
T_eV = 50._dp
T_au = T_ev / Ha2eV
natom = 1
dt = 1e-4_dp
allocate(m(natom))
m = 2._dp * u2au ! Using Argon mass in atomic mass units [u]
velocity_gauge = .true. ! velocity or length gauge?

L = (sum(m) / rho)**(1._dp/3)
L = 10
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
        Ng = 64
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


allocate(ne(Ng_local(1), Ng_local(2), Ng_local(3)))
allocate(neG(Ng_local(1), Ng_local(2), Ng_local(3)))
call allocate_mold(G2, ne)
call allocate_mold(Htot, ne)
call allocate_mold(HtotG, ne)
call allocate_mold(Vloc, ne)
call allocate_mold(Veff, ne)
call allocate_mold(Vxc, ne)
call allocate_mold(exc, ne)
call allocate_mold(Ven0G, ne)
call allocate_mold(Vee, ne)
call allocate_mold(Vee_xc, ne)
call allocate_mold(psiG, neG)
call allocate_mold(psi, neG)
call allocate_mold(tmp, neG)
call allocate_mold(VeeG, neG)
call allocate_mold(VenG, neG)
allocate(X(Ng_local(1), Ng_local(2), Ng_local(3), 3))
allocate(dpsi(Ng_local(1), Ng_local(2), Ng_local(3), 3))
call allocate_mold(G, X)
call allocate_mold(current, X)
if (myid == 0) print *, "Load initial position"
allocate(Xion(3, natom))
allocate(forces(3, natom))
allocate(q(natom))
! For now assume a box, until positions_bcc can accept a vector L(:)
! And radial_potential_fourier
call assert(abs(L(2)-L(1)) < 1e-15_dp)
call assert(abs(L(3)-L(1)) < 1e-15_dp)
!call positions_fcc(Xion, L(1))
!Xion(:,1) = [-0.7_dp+L(1)/2, L(2)/2, L(3)/2]
!Xion(:,2) = [-0.7_dp+L(1)/2, L(2)/2, L(3)/2]
Xion(:,1) = L/2
q = 4
!Xion = 0
!Xion(1,1) = L(1)/2
!Xion = Xion + 1e-4_dp ! Shift the atoms, so that 1/R is defined
call real_space_vectors(L, X, Ng, myxyz)
call reciprocal_space_vectors(L, G, G2, Ng, myxyz)

call ewald_box(L(1), Xion, q, Enn, forces, stress)

! HDH pseudopotential for Hydrogen
rloc = 0.2_dp
C1 = -4.180237_dp
C2 =  0.725075_dp
Zion = 1
! HDH local pseudopotential for Pb
rloc = 0.617500_dp
C1 =   0.753143_dp
C2 =   0
Zion = 4

do k = 1, Ng_local(3)
do j = 1, Ng_local(2)
do i = 1, Ng_local(1)
    w2 = G2(i,j,k)
    if (w2 < tiny(1._dp)) then
        Ven0G(i,j,k) = 0
    else
        Ven0G(i,j,k) = -4*pi*Zion/w2*exp(-w2*rloc**2/2) &
            + sqrt(8*pi**3)*rloc**3*exp(-w2*rloc**2/2)*(C1+C2*(3-w2*rloc**2))
    end if
end do
end do
end do
Ven0G = Ven0G / product(L)


VenG = 0
do i = 1, natom
    VenG = VenG + Ven0G * exp(-i_ * &
        (G(:,:,:,1)*Xion(1,i) + G(:,:,:,2)*Xion(2,i) + G(:,:,:,3)*Xion(3,i)))
end do

!Vloc = 0
!do na = 1, natom
!    do k = 1, Ng_local(3)
!    do j = 1, Ng_local(2)
!    do i = 1, Ng_local(1)
!        r = sqrt(sum((X(i,j,k,:)-Xion(:,na))**2))
!        Vloc(i,j,k) = Vloc(i,j,k) - Zion/r * erf(r/(sqrt(2._dp)*rloc)) &
!            + exp(-1._dp/2*(r/rloc)**2) * (C1 + C2*(r/rloc)**2)
!    end do
!    end do
!    end do
!end do

Vmin = minval(Vloc)

!if (myid == 0) then
!    open(newunit=u, file="a.txt", status="replace")
!    write(u,*) Vloc(:,16,16)
!end if

!call pfourier2real(G2*VenG/(4*pi), Vloc, commy, commz, Ng, nsub)
call pfourier2real(VenG, Vloc, commy, commz, Ng, nsub)
norm = pintegral(comm_all, L, Vloc, Ng)
if (myid == 0) then
    print *, "INT:", norm
end if

!Vloc = Vloc + Vmin - minval(Vloc)

!if (myid == 0) then
!    open(newunit=u, file="a.txt", status="replace")
!    !write(u,*) Vloc(:,16,16)
!    !write(u,*) Ven0G(:,1,1)
!    !write(u,*) Vloc(:,Ng/2,Ng/2)
!    write(u,*) Vloc(:,1,1)
!    close(u)
!end if
!stop "OK"

if (myid == 0) print *, "Solving eigenproblem: DOFs =", product(Ng)

nev = 5
ncv = 100
allocate(eigs(nev), orbitals(Ng_local(1),Ng_local(2),Ng_local(3),nev))
allocate(cutfn(Ng_local(1),Ng_local(2),Ng_local(3)))
cutfn = 1
allocate(occ(4))

occ = [1, 1, 1, 1]
Vee_xc = 0
call mixing_linear(myid, product(Ng_local), Rfunc, 100, 0.7_dp, Vee_xc)


if (myid == 0) print *, "Done"

call mpi_finalize(ierr)

contains

    subroutine Rfunc(x, y, E)
    ! Converge Vee+Vxc only (the other components are constant
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:), E

    ! Schroedinger:
    Veff = Vloc + reshape(x, [Ng_local(1),Ng_local(2),Ng_local(3)])
    call solve_schroedinger(myid, comm_all, commy, commz, Ng, nsub, Veff, &
            L, G2, cutfn, nev, ncv, eigs, orbitals)
    !if (myid == 0) then
    !    print *, "n E"
    !    do i = 1, nev
    !        print *, i, eigs(i)
    !    end do
    !end if

    ! Poisson
    ne = 0
    Ekin = 0
    do i = 1, size(occ)
        ne = ne + occ(i)*orbitals(:,:,:,i)**2
        call preal2fourier(orbitals(:,:,:,i), psiG, commy, commz, Ng, nsub)
        Ekin = Ekin &
            + occ(i) * 1._dp/2 * pintegralG(comm_all, L, G2*abs(psiG)**2)
    end do

    call preal2fourier(ne, neG, commy, commz, Ng, nsub)
    call poisson_kernel(myid, size(neG), neG, G2, VeeG)
    call pfourier2real(VeeG, Vee, commy, commz, Ng, nsub)
    call xc_pz(ne, exc, Vxc)

    Een_loc = pintegral(comm_all, L, Vloc*ne, Ng)
    E_xc = pintegral(comm_all, L, exc*ne, Ng)
    Eee = pintegral(comm_all, L, Vee*ne, Ng) / 2
    Een_core = 4.95047558841102E-02_dp

    Etot = Ekin + Eee + E_xc + Enn + Een_core + Een_loc

    if (myid == 0) then
        do i = 1, nev
            print *, i, eigs(i)
        end do
        print "(a, es22.14)", "Ekin:     ", Ekin
        print "(a, es22.14)", "Eee:      ", Eee
        print "(a, es22.14)", "Exc:      ", E_xc
        print "(a, es22.14)", "Enn:      ", Enn
        print "(a, es22.14)", "Een_core: ", Een_core
        print "(a, es22.14)", "Een_loc:  ", Een_loc
        print "(a, es22.14)", "Een_NL:   ", 0._dp
        print "(a, es22.14)", "Etot:     ", Etot
    end if

    !norm = pintegral(comm_all, L, ne, Ng)
    !if (myid == 0) then
    !    print *, "Density norm:", norm
    !end if

    E = Etot
    y = reshape(Vee + Vxc, [product(Ng_local)]) - x

    end subroutine

end program
