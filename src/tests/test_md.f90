program test_md
use types, only: dp
use constants, only: ang2bohr, K2au, density2gcm3, u2au, s2au
use md, only: velocity_verlet, minimize_energy
use random, only: randn
use utils, only: init_random, stop_error
implicit none

! All variables are in Hartree atomic units

integer, parameter :: N = 864
integer :: i, steps, u
real(dp) :: dt, V(3, N), X(3, N), f(3, N), m(N), L, t, Rcut, eps, sigma, rho
real(dp) :: Temp, Ekin, Epot, Temp_current, t1, t2, t3, t4


steps = 800
dt = 1e-14_dp * s2au
m = 39.948_dp * u2au ! Using Argon mass in atomic mass units [u]
sigma = 3.4_dp * ang2bohr
eps = 120 * K2au
Rcut = 2.25_dp*sigma
rho = 1.374_dp / density2gcm3
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
call random_number(X)
X = L*X

call cpu_time(t1)
call minimize_energy(forces, calc_Epot, X, f, h0=1._dp, max_iter=50)
call cpu_time(t2)

call randn(V)
V = V * sqrt(Temp / spread(m, 1, 3))

open(newunit=u, file="md_results.txt", status="replace", form="unformatted")
t = 0
call cpu_time(t3)
do i = 1, steps
    Ekin = calc_Ekin(V, m)
    Epot = calc_Epot(X)
    Temp_current = 2*Ekin/(3*N)
    print "(i5, ': E=', f10.4, ' a.u.; Epot=', f10.4, ' a.u.; T=',f10.4,' K')",&
        i, Ekin + Epot, Epot, Temp_current / K2au
    write(u) t/s2au, f, V, X, Ekin, Epot, Temp_current / K2au
    call velocity_verlet(dt, m, L, forces, f, V, X)
    t = t + dt
    if (i < 100) then
        Ekin = calc_Ekin(V, m)
        Temp_current = 2*Ekin/(3*N)
        V = V * sqrt(Temp / Temp_current)
    end if
end do
call cpu_time(t4)
close(u)
print "('TIMINGS')"
print "('Minimization:', f10.4, 's')", t2-t1
print "('MD run:      ', f10.4, 's')", t4-t3
print "('MD step:     ', f10.4, 's')", (t4-t3)/steps

contains

    subroutine forces(X, f)
    real(dp), intent(in) :: X(:, :) ! positions
    real(dp), intent(out) :: f(:, :) ! forces
    real(dp) :: r2, d(3), force(3), Xj(3)
    integer :: N, i, j
    N = size(X, 2)
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
