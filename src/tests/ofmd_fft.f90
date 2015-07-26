program ofmd_fft
use, intrinsic :: iso_fortran_env, only: output_unit
use types, only: dp
use constants, only: i_, K2au, density2gcm3, u2au, s2au, Ha2eV
use md, only: velocity_verlet, positions_random, &
                calc_min_distance, positions_fcc
use ewald_sums, only: ewald_box
use random, only: randn
use utils, only: init_random, stop_error, assert, linspace
use ofdft, only: read_pseudo
use ofdft_fft, only: free_energy_min, radial_potential_fourier, &
    reciprocal_space_vectors, real2fourier, fourier2real, logging_info
use interp3d, only: trilinear
use poisson3d_assembly, only: func2quad
implicit none

! All variables are in Hartree atomic units

integer :: N = 4
integer :: i, steps, u
real(dp) :: dt, L, t, rho, scf_eps
real(dp), allocatable :: V(:, :), X(:, :), f(:, :), m(:)
real(dp), allocatable :: R(:), Ven_rad(:), &
    G(:, :, :, :), G2(:, :, :)
real(dp), allocatable :: Ven0G(:, :, :)
complex(dp), allocatable :: VenG(:, :, :), neG(:, :, :)
real(dp), allocatable :: ne(:, :, :), R2(:)
real(dp) :: Temp, Ekin, Temp_current, t3, t4
real(dp) :: Ediff, Z
real(dp) :: Een_correction
real(dp) :: Eee, Een, Ts, Exc, Etot, Enn
real(dp), allocatable :: fnn(:, :), q(:), fen(:, :)
integer :: dynamics, functional, Ng, Nspecies, start, Nmesh

logging_info = .false. ! Turn of the INFO warnings

call read_input("OFMD.input", Temp, rho, Nspecies, N, start, dynamics, &
            functional, Ng, scf_eps, steps, dt)
call read_pseudo("fem/Al.pseudo", R, Ven_rad, Z, Ediff)
allocate(X(3, N), V(3, N), f(3, N), m(N))
allocate(Ven0G(Ng, Ng, Ng), VenG(Ng, Ng, Ng), ne(Ng, Ng, Ng), neG(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng))
allocate(fnn(3, N), q(N), fen(3, N))

q = Z
m = 26.9_dp * u2au ! Using Hydrogen mass in atomic mass units [u]
L = (sum(m) / rho)**(1._dp/3)
print *, "----------------------------------------------------------------"
print *, "Input Summary:"
print *, "N =", N
print *, "Ng =", Ng
print *, "MD steps =", steps
print *, "rho = ",  rho, "a.u. =", rho * density2gcm3, "g/cm^3"
print *, "T =", Temp, "a.u. =", Temp / K2au, "K ="
print *, "  =", Temp * Ha2eV, "eV"
print "(' dt =', f8.2, ' a.u. = ', es10.2, ' s = ', es10.2, ' ps')", dt, &
    dt/s2au, dt/s2au * 1e12_dp
print "(' SCF_eps =', es10.2, ' a.u. = ', es10.2, ' eV')", scf_eps, &
    scf_eps * Ha2eV
print *, "----------------------------------------------------------------"
print *
print *, "Calculated quantities:"
print *, "L =", L, "a.u."
print *

print *, "Converting aperiodic radial Ven to periodic cartesian Ven"
call radial_potential_fourier(R, Ven_rad, L, Z, Ven0G, Een_correction)
print *, "  Done."
Nmesh = 10000
allocate(R2(Nmesh))
R2 = linspace(0._dp, L/2, Nmesh)

call reciprocal_space_vectors(L, G, G2)

! Make it deterministic for now
!call init_random()
!call positions_random(X, L, 2**(1._dp/6)*sigma, 10)
call positions_fcc(X, L)
print *, "Positions:"
do i = 1, N
    print *, i, X(:, i)
end do
print *, "Distances:"
do i = 2, N
    print *, i, calc_min_distance(X(:, :i-1), L, X(:, i))
end do
print *

! Initialize velocities based on Maxwell-Boltzmann distribution
call randn(V)
V = V * sqrt(Temp / spread(m, 1, 3))

Ekin = calc_Ekin(V, m)
Temp_current = 2*Ekin/(3*N)

! The average temperature (i.e. if we average Temp_current for many runs) will
! be Temp. But we want to set it exactly to Temp, so we rescale the velocities.
V = V * sqrt(Temp / Temp_current)

print *, "MD start:"

t = 0
ne = N * Z / L**3

open(newunit=u, file="ofmd_results.txt", status="replace")
call write_results_header(u)
close(u)

t = 0
call cpu_time(t3)
do i = 1, steps
    print *, "Starting MD iteration:", i
    if (i == 1) then
        call forces(X, f)
    else
        call velocity_verlet(dt, m, L, forces, f, V, X)
    end if
    Ekin = calc_Ekin(V, m)
    Temp_current = 2*Ekin/(3*N)

    print *, "Nuclear forces:"
    print *, fnn(:, 1)
    print *, fnn(:, 2)
    print *, fnn(:, 3)
    print *, fnn(:, 4)

    print *, "Electronic forces:"
    print *, fen(:, 1)
    print *, fen(:, 2)
    print *, fen(:, 3)
    print *, fen(:, 4)

    print *, "total forces:"
    print *, f(:, 1)
    print *, f(:, 2)
    print *, f(:, 3)
    print *, f(:, 4)

    open(newunit=u, file="ofmd_results.txt", position="append", status="old")
    call write_results_line(u)
    close(u)

    call write_results_header(output_unit)
    call write_results_line(output_unit)

    t = t + dt
end do
call cpu_time(t4)
close(u)
print "('TIMINGS')"
print "('MD run:      ', f10.4, 's')", t4-t3
print "('MD step:     ', f10.4, 's')", (t4-t3)/steps

contains

    subroutine forces(X, f)
    real(dp), intent(in) :: X(:, :) ! positions
    real(dp), intent(out) :: f(:, :) ! forces
    real(dp) :: stress(6)
    real(dp) :: fac(Ng, Ng, Ng)
    integer :: i
    print *, "-------------------------------------"
    print *, "Calculating forces"
    ! Calculate nuclear forces
    print *, "Calculating nuclear forces"
    call ewald_box(L, X, q, Enn, fnn, stress)

    ! Calculate the electronic forces
    print *, "Calculating VenG"
    VenG = 0
    do i = 1, N
        ! Note: we use minus sign in the forward real -> fourier transform.
        ! Then the following holds:
        !
        !         F[f(x+b)] = F[f(x)]*e^{+i*G*b},
        !
        ! with plus sign in the exponential. Finally, we are expressing
        ! Ven0(x-X) using the above formula (with b=-X) in terms of Ven0(x),
        ! i.e.:
        !
        !         F[Ven0(x-X)] = F[Ven0(x)]*e^{-i*G*X},
        !
        ! with minus sign in the exponential.
        VenG = VenG - Ven0G * exp(i_ * &
            (G(:,:,:,1)*X(1,i) + G(:,:,:,2)*X(2,i) + G(:,:,:,3)*X(3,i)))
    end do
    call assert(abs(VenG(1, 1, 1)) < epsilon(1._dp)) ! The G=0 component

    ! Energy calculation
    print *, "Minimizing free energy"
    call free_energy_min(N*Z, N, L, G2, Temp, VenG, ne, scf_eps, &
            Eee, Een, Ts, Exc, Etot)

    ! Forces calculation
    print *, "ne -> neG"
    call real2fourier(ne, neG)

    Een = Een + Een_correction * real(neG(1, 1, 1), dp) * N
    Etot = Eee + Een + Ts + Exc

    fen = 0
    print *, "Calculating fen"
    do i = 1, N
        ! We have minus sign in the exponential, per the definition of VenG
        ! above (see the comment there).
        fac = L**3*Ven0G*aimag(neG*exp(-i_ * &
            (G(:,:,:,1)*X(1,i) + G(:,:,:,2)*X(2,i) + G(:,:,:,3)*X(3,i))))
        fen(1, i) = sum(G(:,:,:,1)*fac)
        fen(2, i) = sum(G(:,:,:,2)*fac)
        fen(3, i) = sum(G(:,:,:,3)*fac)
    end do

    f = fnn + fen
    print *, "Done calculating forces."
    print *, "-------------------------------------"
    end subroutine

    real(dp) pure function calc_Ekin(V, m) result(Ekin)
    real(dp), intent(in) :: V(:, :) ! velocities
    real(dp), intent(in) :: m(:) ! masses
    Ekin = sum(m*sum(V**2, dim=1))/2
    end function

    subroutine read_input(filename, T, density, Nspecies, N, start, dynamics, &
            functional, Ng, scf_eps, steps, dt)
    ! Reads the input file, returns values in a.u.
    character(len=*), intent(in) :: filename
    real(dp), intent(out) :: T, density, dt, scf_eps
    integer, intent(out) :: Nspecies, N, start, dynamics, functional, Ng, steps
    integer :: u
    real(dp) :: skip
    open(newunit=u, file=filename, status="old")
    read(u, *) T, density, Nspecies, N
    read(u, *) start
    read(u, *) dynamics, functional, Ng, skip, skip, scf_eps
    read(u, *) steps, dt
    close(u)
    ! Convert to atomic units:
    T = T / Ha2eV
    density = density / density2gcm3
    scf_eps = scf_eps / Ha2eV
    end subroutine

    subroutine write_results_header(u)
    integer, intent(in) :: u
    write(u, '(10a17)') "Time [a.u.]", "Fe [eV]", "Unn [eV]", "K [eV]", &
        "F [eV]", "T [eV]", "Ts [eV]", "Een [eV]", "Eee [eV]", "Exc [eV]"
    end subroutine

    subroutine write_results_line(u)
    integer, intent(in) :: u
    write(u, '(10f17.6)') t, Etot*Ha2eV/N, Enn*Ha2eV/N, Ekin*Ha2eV/N, &
        (Etot + Enn + Ekin) * Ha2eV / N, Temp_current * Ha2eV, Ts * Ha2eV / N, &
        Een * Ha2eV / N, Eee * Ha2eV / N, Exc * Ha2eV / N
    end subroutine

end program
