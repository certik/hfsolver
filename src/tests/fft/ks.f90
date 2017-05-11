program ks
use types, only: dp
use constants, only: Ha2eV, density2gcm3, u2au, s2au, K2au, i_, pi, &
    ang2bohr, bohr2ang
use fourier, only: dft, idft, fft, fft_vectorized, fft_pass, fft_pass_inplace, &
        fft_vectorized_inplace, calculate_factors, ifft_pass, fft2_inplace, &
        fft3_inplace, ifft3_inplace
use utils, only: assert, init_random, stop_error, get_int_arg, get_float_arg, &
    allocate_mold, clock, whitechar
use ffte, only: factor
use ofdft, only: read_pseudo
use pofdft_fft, only: pfft3_init, preal2fourier, pfourier2real, &
    real_space_vectors, reciprocal_space_vectors, calculate_myxyz, &
    pintegral, pintegralG, free_energy, free_energy_min, &
    radial_potential_fourier, psum, pmaxval, collate, poisson_kernel
use openmp, only: omp_get_wtime
use mpi2, only: mpi_finalize, MPI_COMM_WORLD, mpi_comm_rank, &
    mpi_comm_size, mpi_init, mpi_comm_split, MPI_INTEGER, &
    mpi_barrier, mpi_bcast, MPI_DOUBLE_PRECISION
use md, only: positions_bcc, positions_fcc
use arpack, only: peig, eig
use pksdft_fft, only: solve_schroedinger
use xc, only: xc_pz
use mixings, only: mixing_linear, mixing_linear_adapt
use ewald_sums, only: ewald_box
use efermi, only: fermi_dirac_smearing
implicit none

complex(dp), dimension(:,:,:), allocatable :: neG, psiG, psi, tmp
real(dp), dimension(:,:,:), allocatable :: G2, Htot, HtotG, Ven0G, ne, Vloc, &
    Veff
real(dp), allocatable :: G(:,:,:,:), X(:,:,:,:), Xion(:,:), q(:), &
    current(:,:,:,:), eigs(:), orbitals(:,:,:,:), eigs_ref(:), occ(:), &
    Vee(:,:,:), Vee_xc(:,:,:), exc(:,:,:), Vxc(:,:,:), forces(:,:), &
    cutfn(:,:,:), cutfn2(:,:,:)
complex(dp), allocatable :: dpsi(:,:,:,:), VeeG(:,:,:), VenG(:,:,:)
real(dp) :: L(3), r, stress(6)
integer :: i, j, k, u
integer :: Ng(3)
integer :: LNPU(3)
integer :: natom, nelec, nband
logical :: velocity_gauge
real(dp) :: T_au, dt, alpha, rho, norm, w2, Vmin, Ekin, Etot, &
    Eee, Een_loc, E_xc, Enn, Een_core, G2cut, G2cut2
real(dp) :: rloc, C1, C2, Zion, Ecut
real(dp), allocatable :: m(:)
integer :: nev, ncv, na
real(dp), parameter :: D(5) = [0.65435_dp, 2.45106_dp, -1.536643785333E-01_dp, &
    1.153664378533E+00_dp, 5.0000_dp]

!  parallel variables
integer :: comm_all, commy, commz, nproc, ierr, nsub(3), Ng_local(3)
integer :: myid ! my ID (MPI rank), starts from 0
integer :: myxyz(3) ! myid, converted to the (x, y, z) box, starts from 0


call mpi_init(ierr)
comm_all  = MPI_COMM_WORLD
call mpi_comm_rank(comm_all, myid, ierr)
call mpi_comm_size(comm_all, nproc, ierr)

if (myid == 0) then
    call load_initial_pos(natom, L, Xion)
end if
call mpi_bcast(natom, 1, MPI_INTEGER, 0, comm_all, ierr)
if (myid /= 0) then
    allocate(Xion(3,natom))
end if
call bcast_float_array(comm_all, size(Xion), Xion)
call bcast_float_array(comm_all, size(L), L)

if (myid == 0) then
    call read_input(nproc, Ng, nsub, T_au, dt, Ecut, nband)
end if
call mpi_bcast(Ng, size(Ng), MPI_INTEGER, 0, comm_all, ierr)
call mpi_bcast(nband, 1, MPI_INTEGER, 0, comm_all, ierr)
call mpi_bcast(nsub, size(nsub), MPI_INTEGER, 0, comm_all, ierr)
call mpi_bcast(T_au, 1, MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
call mpi_bcast(dt, 1, MPI_DOUBLE_PRECISION, 0, comm_all, ierr)
call mpi_bcast(Ecut, 1, MPI_DOUBLE_PRECISION, 0, comm_all, ierr)

G2cut = 2*Ecut
G2cut2 = 4*G2cut
if (any(Ng == -1)) then
    Ng(1) = int(2*sqrt(G2cut2) * L(1) / (2*pi) + 1)
    call adjust_Ng(Ng(1))
    Ng(2:3) = Ng(1)
end if

Ng_local = Ng / nsub


allocate(m(natom))
nelec = 4*natom
m = 2._dp * u2au ! Using D mass in atomic mass units [u]
velocity_gauge = .true. ! velocity or length gauge?

!L = (sum(m) / rho)**(1._dp/3)
rho = sum(m) / product(L)

allocate(q(natom))
q = 1

allocate(occ(nband))

call save_abinit(L, Xion, T_au, dt, Ecut, m, q)

if (myid == 0) then
    print *, "Input:"
    print *, "N =", natom
    print *, "L =", L, "a.u. =", L * bohr2ang, "Angst"
    print *, "T =", T_au / K2au, "K =", T_au, "a.u. =", T_au * Ha2eV, "eV"
    print "('dt =', es10.2, ' s = ', es10.2, ' fs = ', es10.2, ' a.u.')", &
        dt/s2au, dt/s2au * 1e15_dp, dt
    print *, "Ecut =", Ecut, "a.u. =", Ecut * Ha2eV, "eV"
    print *, "nband =", nband
    print *
    print *, "Calculated quantities:"
    print *, "rho = ", rho * density2gcm3, "g/cc = ", rho, "a.u."
    print *
    print *, "nproc:   ", nproc
    print *, "nsub:    ", nsub
    print *, "Ng:      ", Ng
    print *, "Ng_local:", Ng_local

    if (product(nsub) /= nproc) then
        call stop_error("nproc must be equal to the number of subdomains")
    end if
end if

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
allocate(forces(3, natom))
! For now assume a box, until positions_bcc can accept a vector L(:)
! And radial_potential_fourier
call assert(abs(L(2)-L(1)) < 1e-15_dp)
call assert(abs(L(3)-L(1)) < 1e-15_dp)
!call positions_fcc(Xion, L(1))
!Xion(:,1) = [-0.7_dp+L(1)/2, L(2)/2, L(3)/2]
!Xion(:,2) = [-0.7_dp+L(1)/2, L(2)/2, L(3)/2]
!Xion(:,1) = L/2
!Xion = 0
!Xion(1,1) = L(1)/2
!Xion = Xion + 1e-4_dp ! Shift the atoms, so that 1/R is defined
call real_space_vectors(L, X, Ng, myxyz)
call reciprocal_space_vectors(L, G, G2, Ng, myxyz)

call ewald_box(L(1), Xion, q, Enn, forces, stress)

! HDH pseudopotential for Hydrogen
!rloc = 0.2_dp
!C1 = -4.180237_dp
!C2 =  0.725075_dp
!Zion = 1
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

!print *, "Ecut:", Ecut
!print *, "G2cut:", G2cut / (2*pi)**2
!print *, "boxsq:", maxval(G2(:,1,1)) / (2*pi)**2
!print *, "boxcut:", sqrt(maxval(G2(:,1,1)) / G2cut)
!print *, "gsqcut:", G2cut2 / (2*pi)**2

allocate(cutfn(Ng_local(1),Ng_local(2),Ng_local(3)))
allocate(cutfn2(Ng_local(1),Ng_local(2),Ng_local(3)))
where (G2 > G2cut)
    cutfn = 0
elsewhere
    !cutfn = (1-G2/G2cut)**12
    cutfn = 1
end where
!print *, maxval(G2), G2cut, 2.04204_dp**2 * G2cut
!stop "OK"
where (G2 > G2cut2)
    cutfn2 = 0
elsewhere
    cutfn2 = 1
end where

!call pfourier2real(G2*VenG/(4*pi), Vloc, commy, commz, Ng, nsub)
VenG = VenG * cutfn2
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

nev = nband
ncv = max(4*nev, 100)
allocate(eigs(nev), orbitals(Ng_local(1),Ng_local(2),Ng_local(3),nev))
Vee_xc = 0
call mixing_linear(myid, product(Ng_local), Rfunc, 100, 0.7_dp, Vee_xc)


if (myid == 0) print *, "Done"

call mpi_finalize(ierr)

contains

    subroutine Rfunc(x, y, E)
    ! Converge Vee+Vxc only (the other components are constant
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:), E
    real(dp) :: mu, sigma

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

    sigma = T_au
    call fermi_dirac_smearing(eigs, sigma, real(nelec, dp), mu, occ)

    ! Poisson
    ne = 0
    Ekin = 0
    do i = 1, nband
        ne = ne + occ(i)*orbitals(:,:,:,i)**2
        call preal2fourier(orbitals(:,:,:,i), psiG, commy, commz, Ng, nsub)
        Ekin = Ekin &
            + occ(i) * 1._dp/2 * pintegralG(comm_all, L, G2*abs(psiG)**2)
    end do

    call preal2fourier(ne, neG, commy, commz, Ng, nsub)
    call poisson_kernel(myid, size(neG), neG, G2, VeeG)
    VeeG = VeeG * cutfn2
    call pfourier2real(VeeG, Vee, commy, commz, Ng, nsub)
    call xc_pz(ne, exc, Vxc)

    Een_loc = pintegral(comm_all, L, Vloc*ne, Ng)
    E_xc = pintegral(comm_all, L, exc*ne, Ng)
    Eee = pintegral(comm_all, L, Vee*ne, Ng) / 2
    Een_core = 4.95047558841102E-02_dp

    Etot = Ekin + Eee + E_xc + Enn + Een_core + Een_loc

    if (myid == 0) then
        !eigs = eigs * Ha2eV
        print *, "E_fermi =", mu
        do i = 1, nband
            print *, i, eigs(i), occ(i)
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

    subroutine read_input(nproc, Ng, nsub, T, dt, Ecut, nband)
    integer, intent(in) :: nproc
    integer, intent(out) :: Ng(3), nsub(3), nband
    real(dp), intent(out) :: T  ! in a.u.
    real(dp), intent(out) :: dt ! in a.u.
    real(dp), intent(out) :: ecut ! in a.u.
    integer :: LNPU(3)
    real(dp) :: T_eV, ecut_eV
    namelist /domain/ Ng, nsub, T_eV, dt, ecut, ecut_eV, nband
    integer :: u
    Ng = -1
    T_eV = -1
    dt = -1
    ecut = -1
    ecut_eV = -1
    nband = -1
    open(newunit=u, file="input", status="old")
    read(u,nml=domain)
    close(u)
    if (nsub(1) == -1) then
        call factor(nproc, LNPU)
        nsub(3) = (2**(LNPU(1)/2))*(3**(LNPU(2)/2))*(5**(LNPU(3)/2))
        nsub(2) = nproc / nsub(3)
        nsub(1) = 1
    end if
    if (Ng(2) == -1 .and. Ng(3) == -1) Ng(2:3) = Ng(1)
    if (T_eV < 0) call stop_error("T is not specified")
    if (dt < 0) call stop_error("dt is not specified")
    if (ecut < 0) then
        if (ecut_eV < 0) call stop_error("ecut nor ecut_eV not specified")
        ecut = ecut_eV / Ha2eV ! Convert from eV to a.u.
    end if
    if (nband < 0) call stop_error("Must specify nband")

    T = T_eV / Ha2eV  ! Convert from eV to a.u.
    endsubroutine

    subroutine load_initial_pos(natom, L, Xion)
    integer, intent(out) :: natom
    real(dp), intent(out) :: L(:)
    real(dp), intent(out), allocatable :: Xion(:,:)
    real(dp) :: A, t(3,3)
    integer :: u, ios, natom_types
    character :: c
    integer, allocatable :: ncounts(:)
    open(newunit=u, file="POSCAR", status="old")
    read(u,*) ! 1 line description
    read(u,*) A
    read(u,*) t(:, 1)
    read(u,*) t(:, 2)
    read(u,*) t(:, 3)

    t = t*A*ang2bohr
    L(1) = t(1,1)
    L(2) = t(2,2)
    L(3) = t(3,3)
    call assert(abs(t(2, 1)) < tiny(1._dp))
    call assert(abs(t(3, 1)) < tiny(1._dp))
    call assert(abs(t(1, 2)) < tiny(1._dp))
    call assert(abs(t(3, 2)) < tiny(1._dp))
    call assert(abs(t(1, 3)) < tiny(1._dp))
    call assert(abs(t(2, 3)) < tiny(1._dp))

    natom_types = number_of_columns(u)
    allocate(ncounts(natom_types))
    read(u,*) ! Atom types
    read(u,*) ncounts ! Atom numbers per type
    natom = sum(ncounts)
    allocate(Xion(3,natom))
    read(u,*) c ! Selective dynamics
    call assert(c == "S")
    read(u,*) c ! Direct
    call assert(c == "D")
    do i = 1, natom
        read(u,*) Xion(:,i)
        Xion(:,i) = Xion(:,i)*L
    end do
    close(u)
    end subroutine

    integer function number_of_columns(u) result(ncol)
    ! Determines the number of columns on a line (does not affect file
    ! position).
    integer, intent(in) :: u ! file handle
    integer :: ios
    character :: c
    logical :: lastwhite
    ncol = 0
    lastwhite = .true.
    do
       read(u, '(a)', advance='no', iostat=ios) c
       if (ios /= 0) exit
       if (lastwhite .and. .not. whitechar(c)) ncol = ncol + 1
       lastwhite = whitechar(c)
    end do
    backspace(u)
    end function

    subroutine bcast_float_array(comm, n, a)
    integer, intent(in) :: comm, n
    real(dp), intent(inout) :: a(n)
    integer :: ierr
    call mpi_bcast(a, n, MPI_DOUBLE_PRECISION, 0, comm, ierr)
    end subroutine

    subroutine adjust_Ng(Ng)
    integer, intent(inout) :: Ng
    integer :: LNPU(3)
    do while (.true.)
        call factor(Ng, LNPU)
        if ((2**LNPU(1))*(3**LNPU(2))*(5**LNPU(3)) == Ng) return
        Ng = Ng + 1
    end do
    end subroutine

    subroutine save_abinit(L, Xion, T_au, dt, Ecut, m, q)
    real(dp), intent(in) :: L(:), Xion(:,:), T_au, dt, Ecut, m(:), q(:)
    integer :: u, natom, i, typat(size(Xion,2))
    natom = size(Xion,2)
    typat = 1
    open(newunit=u, file="abinit.in", status="replace")
    write(u,*) "acell", L
    write(u,*) "ntypat", 1
    write(u,*) "znucl", 1 ! Hydrogen atom type
    write(u,*) "natom", natom
    write(u,*) "typat"
    do i = 1, natom
        write(u,*) typat(i)
    end do
    write(u,*) "xcart"
    do i = 1, natom
        write(u,*) Xion(:,i)
    end do
    write(u,*) "ecut", Ecut
    write(u,*) "kptopt", 0
    write(u,*) "nkpt", 1
    write(u,*) "nstep", 50
    write(u,*) "toldfe", 1e-12
    write(u,*) "diemac", 2._dp
    write(u,*) "occopt", 3 ! Fermi-Dirac smearing
    write(u,*) "tsmear", T_au
    write(u,*) "nband", size(occ)
    write(u,*) "ixc", 2
    write(u,*) "istwfk", 1
    write(u,*) "nsym", 1
    close(u)
    end subroutine

end program
