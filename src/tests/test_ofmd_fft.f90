program test_ofmd_fft
use types, only: dp
use constants, only: ang2bohr, K2au, density2gcm3, u2au, s2au
use md, only: velocity_verlet, minimize_energy, random_positions, &
                calc_min_distance
use ewald_sums, only: ewald
use random, only: randn
use utils, only: init_random, stop_error
implicit none

! All variables are in Hartree atomic units

integer, parameter :: N = 4
integer :: i, steps, u
real(dp) :: dt, V(3, N), X(3, N), f(3, N), m(N), L, t, Rcut, eps, sigma, rho
real(dp) :: Temp, Ekin, Epot, Temp_current, t3, t4


steps = 800
dt = 1e-14_dp * s2au
m = 39.948_dp * u2au ! Using Argon mass in atomic mass units [u]
sigma = 3.4_dp * ang2bohr
eps = 120 * K2au
Rcut = 2.25_dp*sigma
!rho = 1.374_dp / density2gcm3
rho = 0.4_dp / density2gcm3
Temp = 94.4_dp * K2au

L = (sum(m) / rho)**(1._dp/3)
rho = sum(m) / L**3
print *, "Input:"
print *, "N =", N
print *, "rho = ", rho * density2gcm3, "g/cm^3 = ", rho, "a.u."
print *, "T =", Temp / K2au, "K =", Temp, "a.u."
print "('dt =', es10.2, ' s = ', es10.2, ' ps = ', f8.2, ' a.u.')", dt/s2au, &
    dt/s2au * 1e12_dp, dt
print *, "eps =", eps / K2au, "K/kB =", eps, "a.u.";
print *, "sigma =", sigma / ang2bohr, " A =", sigma, "a.u."
print *, "Rcut =", Rcut / sigma, " sigma =", Rcut, "a.u."
print *
print *, "Calculated quantities:"
print *, "L =", L/sigma, "sigma = ", L, "a.u."
print *

call init_random()
call random_positions(X, L, 2**(1._dp/6)*sigma, 10)
print *, "Distances:"
do i = 2, N
    print *, i, calc_min_distance(X(:, :i-1), L, X(:, i))
end do
print *

call randn(V)
V = V * sqrt(Temp / spread(m, 1, 3))

print *, "MD start:"

call forces(X, f)

open(newunit=u, file="ofmd_results.txt", status="replace", form="unformatted")
t = 0
call cpu_time(t3)
do i = 1, steps
    Ekin = calc_Ekin(V, m)
    Epot = calc_Epot(X)
    Temp_current = 2*Ekin/(3*N)
    print "(i5, ': E=', f10.4, ' a.u.; Epot=', f10.4, ' a.u.; T=',f10.4,' K')",&
        i, Ekin + Epot, Epot, Temp_current / K2au
    write(u) t/s2au, f, V, X, Ekin, Epot, Temp_current / K2au
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
    real(dp) :: r2, d(3), force(3), Xj(3)
    integer :: N, i, j
    integer :: ntypat
    real(dp) :: ucvol
    integer, allocatable :: typat(:)
    real(dp) :: gmet(3, 3), rmet(3, 3)
    real(dp) :: eew
    real(dp), allocatable :: xred(:, :), zion(:), grewtn(:, :)
    N = size(X, 2)
    ntypat = 1
    ! TODO: this can be done in the main program
    allocate(xred(3, N), zion(ntypat), grewtn(3, N), typat(N))
    xred = X / L
    typat = 1
    zion = [1._dp]

    rmet = 0
    rmet(1, 1) = L**2
    rmet(2, 2) = L**2
    rmet(3, 3) = L**2

    ! gmet = inv(rmet)
    gmet = 0
    gmet(1, 1) = 1/L**2
    gmet(2, 2) = 1/L**2
    gmet(3, 3) = 1/L**2

    ! ucvol = sqrt(det(rmet))
    ucvol = L**3

    call ewald(eew,gmet,grewtn,N,ntypat,rmet,typat,ucvol,xred,zion)

    f = 0
    do i = 1, N
        do j = 1, i-1
            Xj = X(:, j)-X(:, i)+[L/2, L/2, L/2]
            Xj = Xj - L*floor(Xj/L)
            d = [L/2, L/2, L/2] - Xj
            r2 = sum(d**2)
            if (r2 <= Rcut**2) then
                force = 24*eps*d/r2 * (2*sigma**12/r2**6 - sigma**6/r2**3)
                f(:, i) = f(:, i) + force
                f(:, j) = f(:, j) - force
            end if
        end do
    end do
    end subroutine

    real(dp) function calc_Epot(X) result(E)
    real(dp), intent(in) :: X(:, :) ! positions
    real(dp) :: r2, d(3), Xj(3)
    integer :: N, i, j
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

end program
