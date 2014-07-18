program test_ofmd_fft
use types, only: dp
use constants, only: i_, K2au, density2gcm3, u2au, s2au, Ha2eV
use md, only: velocity_verlet, minimize_energy, positions_random, &
                calc_min_distance, positions_fcc
use ewald_sums, only: ewald_box
use random, only: randn
use utils, only: init_random, stop_error, assert
use ofdft, only: read_pseudo
use ofdft_fft, only: free_energy_min, radial_density_fourier, &
    reciprocal_space_vectors, real2fourier
implicit none

! All variables are in Hartree atomic units

integer :: N = 4
integer :: i, steps, u
real(dp) :: dt, L, t, Rcut, eps, sigma, rho
real(dp), allocatable :: V(:, :), X(:, :), f(:, :), m(:)
real(dp), allocatable :: R(:), Ven_rad(:), &
    G(:, :, :, :), G2(:, :, :)
real(dp), allocatable :: Ven0G(:, :, :)
complex(dp), allocatable :: VenG(:, :, :), neG(:, :, :)
real(dp), allocatable :: ne(:, :, :)
real(dp) :: Temp, Ekin, Epot, Temp_current, t3, t4
real(dp) :: Ediff, Z
integer :: dynamics, functional, Ng, Nspecies, start

call read_input("OFMD.input", Temp, rho, Nspecies, N, start, dynamics, &
            functional, Ng, steps, dt)
call read_pseudo("fem/H.pseudo.gaussian", R, Ven_rad, Z, Ediff)
allocate(X(3, N), V(3, N), f(3, N), m(N))
allocate(Ven0G(Ng, Ng, Ng), VenG(Ng, Ng, Ng), ne(Ng, Ng, Ng), neG(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng))

m = 1._dp * u2au ! Using Hydrogen mass in atomic mass units [u]
L = (sum(m) / rho)**(1._dp/3)
print *, "Input:"
print *, "N =", N
print *, "Ng =", Ng
print *, "MD steps =", steps
print *, "rho = ",  rho, "a.u. =", rho * density2gcm3, "g/cm^3"
print *, "T =", Temp, "a.u. =", Temp / K2au, "K =", Temp * Ha2eV, "eV"
print "('dt =', f8.2, ' a.u. = ', es10.2, ' s = ', es10.2, ' ps')", dt, &
    dt/s2au, dt/s2au * 1e12_dp
print *
print *, "Calculated quantities:"
print *, "L =", L, "a.u."
print *

call radial_density_fourier(R, Ven_rad, L, Z, Ven0G)

call reciprocal_space_vectors(L, G, G2)

call init_random()
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

call randn(V)
V = V * sqrt(Temp / spread(m, 1, 3))

print *, "MD start:"

t = 0
open(newunit=u, file="ofmd_results.txt", status="replace")
ne=1._dp / L**3
call forces(X, f)

t = 0
call cpu_time(t3)
do i = 1, steps
    print *, "MD step:", i, Temp, Temp_current
    Ekin = calc_Ekin(V, m)
    Epot = calc_Epot(X)
    Temp_current = 2*Ekin/(3*N)
    print "(i5, ': E=', f10.4, ' a.u.; Epot=', f10.4, ' a.u.; T=',f10.4,' K')",&
        i, Ekin + Epot, Epot, Temp_current / K2au
    call velocity_verlet(dt, m, forces, f, V, X)
    if (any(abs(X/L) > 2)) then
        print *, "max n = X/L =", maxval(abs(X/L))
        call stop_error("X is out of range after Verlet: abs(X/L) > 2")
    end if
    ! Periodically shift particles to the [0, L]^3 box
    X = X - L*floor(X/L)
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
    integer :: N
    real(dp) :: stress(6)
    real(dp) :: E_ewald
    real(dp), allocatable :: fewald(:, :), q(:), fen(:, :)
    real(dp) :: fac(Ng, Ng, Ng)
    real(dp) :: Eee, Een, Ts, Exc, Etot
    integer :: i
    N = size(X, 2)
    ! TODO: this can be done in the main program
    allocate(fewald(3, N), q(N), fen(3, N))
    q = 1
    ! Calculate nuclear forces
    call ewald_box(L, X, q, E_ewald, fewald, stress)

    ! Calculate the electronic forces
    VenG = 0
    do i = 1, N
        VenG = VenG - Ven0G * exp(i_ * &
            (G(:,:,:,1)*X(1,i) + G(:,:,:,2)*X(2,i) + G(:,:,:,3)*X(3,i)))
    end do
    call assert(abs(VenG(1, 1, 1)) < epsilon(1._dp)) ! The G=0 component
    call free_energy_min(N, L, G2, Temp, VenG, ne, Eee, Een, Ts, Exc, Etot)
    print *, "Ng =", Ng
    print *, "Rcut =", Rcut
    print *, "T_au =", Temp
    print *, "Summary of energies [a.u.]:"
    print "('    Ts   = ', f14.8)", Ts
    print "('    Een  = ', f14.8)", Een
    print "('    Eee  = ', f14.8)", Eee
    print "('    Exc  = ', f14.8)", Exc
    print *, "   ---------------------"
    print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV

    write(u, *) t, Etot*Ha2eV/N, E_ewald*Ha2eV/N, Ekin*Ha2eV/N, Temp_current / K2au

    print *, "EWALD", E_ewald
    print *, fewald(:, 1)
    print *, fewald(:, 2)
    print *, fewald(:, 3)
    print *, fewald(:, 4)
    print *, "stress"
    print *, stress

    call real2fourier(ne, neG)

    fen = 0
    do i = 1, N
        fac = L**3*Ven0G*aimag(neG*exp(-i_ * &
            (G(:,:,:,1)*X(1,i) + G(:,:,:,2)*X(2,i) + G(:,:,:,3)*X(3,i))))
        fen(1, i) = sum(G(:,:,:,1)*fac)
        fen(2, i) = sum(G(:,:,:,2)*fac)
        fen(3, i) = sum(G(:,:,:,3)*fac)
    end do

    print *, "forces:"
    print *, fen(:, 1)
    print *, fen(:, 2)
    print *, fen(:, 3)
    print *, fen(:, 4)


    f = fewald + fen

    print *, "total forces:"
    print *, f(:, 1)
    print *, f(:, 2)
    print *, f(:, 3)
    print *, f(:, 4)
    end subroutine

    real(dp) function calc_Epot(X) result(E)
    real(dp), intent(in) :: X(:, :) ! positions
    real(dp) :: r2, d(3), Xj(3)
    integer :: N, i, j
    ! TODO: this needs to get calculated as part of the force calculation
    N = size(X, 2)
    E = 0
    do i = 1, N
        do j = 1, i-1
            Xj = X(:, j)-X(:, i)+[L/2, L/2, L/2]
            Xj = Xj - L*floor(Xj/L)
            d = [L/2, L/2, L/2] - Xj
            r2 = sum(d**2)
            if (r2 <= Rcut**2) E = E + 4*eps *(sigma**12/r2**6 - sigma**6/r2**3)
        end do
    end do
    end function

    real(dp) pure function calc_Ekin(V, m) result(Ekin)
    real(dp), intent(in) :: V(:, :) ! velocities
    real(dp), intent(in) :: m(:) ! masses
    Ekin = sum(m*sum(V**2, dim=1))/2
    end function

    subroutine read_input(filename, T, density, Nspecies, N, start, dynamics, &
            functional, Ng, steps, dt)
    ! Reads the input file, returns values in a.u.
    character(len=*), intent(in) :: filename
    real(dp), intent(out) :: T, density, dt
    integer, intent(out) :: Nspecies, N, start, dynamics, functional, Ng, steps
    integer :: u
    open(newunit=u, file=filename, status="old")
    read(u, *) T, density, Nspecies, N
    read(u, *) start
    read(u, *) dynamics, functional, Ng
    read(u, *) steps, dt
    close(u)
    ! Convert to atomic units:
    T = T / Ha2eV
    density = density / density2gcm3
    end subroutine

end program
