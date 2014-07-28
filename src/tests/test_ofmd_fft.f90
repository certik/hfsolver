program test_ofmd_fft
use types, only: dp
use constants, only: i_, K2au, density2gcm3, u2au, s2au, Ha2eV, pi
use md, only: velocity_verlet, minimize_energy, positions_random, &
                calc_min_distance, positions_fcc
use ewald_sums, only: ewald_box
use random, only: randn
use utils, only: init_random, stop_error, assert, linspace
!use feutils, only: quad_lobatto
use feutils, only: quad_gauss
use ofdft, only: read_pseudo
use ofdft_fft, only: free_energy_min, radial_potential_fourier, &
    reciprocal_space_vectors, real2fourier, fourier2real
use ofdft_fe, only: radial_density_fourier, initialize_fe, &
    free_energy_min_low_level, fe_data
use interp3d, only: trilinear
use poisson3d_assembly, only: func2quad
use fe_mesh, only: quad2fe_3d, fe_eval_xyz
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
real(dp), allocatable :: nen(:, :, :)
real(dp), allocatable :: ne(:, :, :), R2(:), nen0(:), fullsol(:)
real(dp) :: Temp, Ekin, Epot, Temp_current, t3, t4
real(dp) :: Ediff, Z
integer :: dynamics, functional, Ng, Nspecies, start, Nmesh
integer :: Nx, Ny, Nz, p, Nq, quad_type
type(fe_data) :: fed
real(dp), allocatable, dimension(:, :, :, :) :: nenq_pos, nq_pos


call read_input("OFMD.input", Temp, rho, Nspecies, N, start, dynamics, &
            functional, Ng, steps, dt)
call read_pseudo("fem/H.pseudo.gaussian2", R, Ven_rad, Z, Ediff)
allocate(X(3, N), V(3, N), f(3, N), m(N))
allocate(Ven0G(Ng, Ng, Ng), VenG(Ng, Ng, Ng), ne(Ng, Ng, Ng), neG(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng))
allocate(nen(Ng, Ng, Ng))

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

call radial_potential_fourier(R, Ven_rad, L, Z, Ven0G)
Nmesh = 10000
allocate(R2(Nmesh), nen0(Nmesh))
R2 = linspace(0._dp, L/2, Nmesh)
call radial_density_fourier(R, Ven_rad, L, Z, Ng, R2, nen0)
open(newunit=u, file="H.pseudo.density", status="replace")
write(u, "(a)") "# Density. The lines are: r, nen(r)"
write(u, *) R2
write(u, *) nen0
close(u)

Nx = 3
Ny = 3
Nz = 3
p = 6
Nq = 30
quad_type = quad_gauss
call initialize_fe(L, Nx, Ny, Nz, p, Nq, quad_type, fed)
allocate(nenq_pos(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(nq_pos(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(fullsol(maxval(fed%in)))
nq_pos = func2quad(fed%nodes, fed%elems, fed%xiq, ne_fn)

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
    real(dp) :: Eee_fe, Een_fe, Ts_fe, Exc_fe, Etot_fe
    integer :: i, j, k
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

    ! Energy calculation
    ! old eps: 3.6749308286427368e-5_dp
    call free_energy_min(N, L, G2, Temp, VenG, ne, 1e-9_dp, &
            Eee, Een, Ts, Exc, Etot)
    print *, "Ng =", Ng
    print *, "Rcut =", Rcut
    print *, "T_au =", Temp
    print *, "Summary of FFT energies [a.u.]:"
    print "('    Ts   = ', f14.8)", Ts
    print "('    Een  = ', f14.8)", Een
    print "('    Eee  = ', f14.8)", Eee
    print "('    Exc  = ', f14.8)", Exc
    print *, "   ---------------------"
    print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV


    ! TODO: The nen should rather be calculated directly using nen0
    call fourier2real(VenG*G2/(4*pi), nen)
    nenq_pos = func2quad(fed%nodes, fed%elems, fed%xiq, nen_fn)

    call free_energy_min_low_level(real(N, dp), Temp, nenq_pos, nq_pos, &
        1e-9_dp, fed, Eee_fe, Een_fe, Ts_fe, Exc_fe)
    Etot_fe = Eee_fe + Een_fe + Ts_fe + Exc_fe
    print *, "Summary of FE energies [a.u.]:"
    print "('    Ts   = ', f14.8)", Ts_fe
    print "('    Een  = ', f14.8)", Een_fe
    print "('    Eee  = ', f14.8)", Eee_fe
    print "('    Exc  = ', f14.8)", Exc_fe
    print *, "   ---------------------"
    print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot_fe, &
        Etot_fe*Ha2eV
    write(u, *) t, Etot*Ha2eV/N, E_ewald*Ha2eV/N, Ekin*Ha2eV/N, &
        Temp_current / K2au, Etot_fe*Ha2eV/N

    print *, "EWALD", E_ewald
    print *, fewald(:, 1)
    print *, fewald(:, 2)
    print *, fewald(:, 3)
    print *, fewald(:, 4)
    print *, "stress"
    print *, stress

    ! Forces calculation
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

    ! Calculate forces using FE density
    print *, "FULLSOL"
    call quad2fe_3d(fed%Ne, fed%Nb, fed%p, fed%jac_det, fed%wtq3, &
            fed%Sp, fed%Sj, fed%Sx, fed%phi_v, fed%in, fed%ib, &
            nq_pos, fullsol)
    print *, "ne (FFT) ="
    print *, ne(:3, :3, :3)
    print *, "eval uniform grid"
    do i = 1, Ng
    do j = 1, Ng
    do k = 1, Ng
        ne(i, j, k) = fe_eval_xyz(fed%xin, fed%nodes, fed%elems, fed%in, &
            fullsol, [L, L, L]/(Ng+1) * [i, j, k])
    end do
    end do
    end do
    print *, "Done"
    print *, "ne (FE) ="
    print *, ne(:3, :3, :3)

    call real2fourier(ne, neG)

    fen = 0
    do i = 1, N
        fac = L**3*Ven0G*aimag(neG*exp(-i_ * &
            (G(:,:,:,1)*X(1,i) + G(:,:,:,2)*X(2,i) + G(:,:,:,3)*X(3,i))))
        fen(1, i) = sum(G(:,:,:,1)*fac)
        fen(2, i) = sum(G(:,:,:,2)*fac)
        fen(3, i) = sum(G(:,:,:,3)*fac)
    end do

    print *, "forces FE:"
    print *, fen(:, 1)
    print *, fen(:, 2)
    print *, fen(:, 3)
    print *, fen(:, 4)
    stop "OK"


    f = fewald + fen

    print *, "total forces:"
    print *, f(:, 1)
    print *, f(:, 2)
    print *, f(:, 3)
    print *, f(:, 4)
    end subroutine

    real(dp) function nen_fn(x, y, z) result(n)
    real(dp), intent(in) :: x, y, z
    n = trilinear([x, y, z], [0, 0, 0]*1._dp, [L, L, L], nen)
    end function

    real(dp) function ne_fn(x, y, z) result(n)
    real(dp), intent(in) :: x, y, z
    n = x+y+z ! Silence compiler warning
    n = 1
    end function

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
