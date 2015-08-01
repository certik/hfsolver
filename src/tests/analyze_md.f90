program analyze_md
use types, only: dp
use constants, only: ang2bohr, K2au, density2gcm3, u2au, s2au
use md, only: unfold_positions, pair_correlation_function
use random, only: randn
use utils, only: init_random
use optimize, only: linregress
implicit none

! All variables are in Hartree atomic units

integer, parameter :: N = 864, Daverages=64
integer :: i, j, steps, u, N0, Nsteps, N0list(Daverages)
real(dp) :: dt, f(3, N), m(N), Rcut, eps, sigma, rho, &
    N0list_r(Daverages)
real(dp) :: Temp, Ekin, Epot, s, L
real(dp), allocatable :: X(:, :, :), V(:, :, :), xx(:), yy(:), velcorr(:), &
    Xu(:, :, :), com(:, :), Temp_current(:), t(:), R(:), gr(:)
real(dp) :: slope, intercept, r_, stderr_slope, stderr_intercept


steps = 800
allocate(X(3, N, steps), Xu(3, N, steps), t(steps))
allocate(V(3, N, steps), com(3, steps), Temp_current(steps))
dt = 1e-14_dp * s2au
m = 39.948_dp * u2au ! Using Argon mass in atomic mass units [u]
sigma = 3.4_dp * ang2bohr
eps = 120 * K2au
Rcut = 2.25_dp*sigma
rho = 1.374_dp / density2gcm3
Temp = 94.4_dp * K2au
L = 65.537399673372406_dp

print *, "Loading data..."
open(newunit=u, file="md_results.txt", status="old", form="unformatted")
do i = 1, steps
    read(u) t(i), f, V(:, :, i), X(:, :, i), Ekin, Epot, Temp_current(i)
end do
close(u)
print *, "Done."

print *, "Unfolding positions..."
call unfold_positions(L, X, Xu)
print *, "Done."
do i = 1, steps
    com(:, i) = sum(spread(m, 1, 3) * Xu(:, :, i), dim=2) / sum(m)
end do
do i = 2, steps
    com(:, i) = com(:, i) - com(:, 1)
end do
com(:, 1) = 0
do i = 1, N
    Xu(:, i, :) = Xu(:, i, :) - com(:, :)
end do

N0 = 200
print *, "Tavg =", sum(Temp_current(N0:)) / (steps-N0+1)

Nsteps = 300
call init_random()
call random_number(N0list_r)
N0list = int(200 + 300*N0list_r)
allocate(xx(Nsteps), yy(Nsteps), velcorr(Nsteps))
do i = 1, Nsteps
    s = 0
    do j = 1, size(N0list)
        N0 = N0list(j)
        s = s + sum((Xu(:, :, i+N0) - Xu(:, :, N0))**2) / N
    end do
    s = s / size(N0list)
    xx(i) = t(i)
    yy(i) = s

    do j = 1, size(N0list)
        N0 = N0list(j)
        s = sum(V(:, :, i+N0) * V(:, :, N0)) / N
    end do
    s = s / size(N0list)
    velcorr(i) = s
end do
call linregress(xx(51:), yy(51:)*(0.529177249_dp*1e-8)**2, &
    slope, intercept, r_, stderr_slope, stderr_intercept)
slope = slope / 6
stderr_slope = stderr_slope / 6
print *, "Diffusion:"
call print_valerr(slope, stderr_slope)
print *, "intercept:"
call print_valerr(intercept, stderr_intercept)
open(newunit=u, file="D.txt", status="replace")
do i = 1, size(xx)
    write(u, *) xx(i), yy(i), velcorr(i), com(:, 200+i)
end do
close(u)

print *, "Temperature"
call linregress(t(400:)*1e12_dp, Temp_current(400:), &
    slope, intercept, r_, stderr_slope, stderr_intercept)
print *, "Slope:"
call print_valerr(slope, stderr_slope)
print *, "intercept:"
call print_valerr(intercept, stderr_intercept)

print *, "Pair correlation function"
allocate(R(200))
allocate(gr(size(R)))
call pair_correlation_function(L, X, R, gr)

open(newunit=u, file="g_r.txt", status="replace")
do i = 1, size(R)
    write (u, *)  R(i), gr(i)
end do
close(u)

contains

subroutine print_valerr(val, err)
real(dp), intent(in) :: val, err
integer :: e
e = getexp(val)
print "('(', f0.6, ' +/- ', f0.6, ') 1e', i0)", val*10._dp**(-e), &
    err*10._dp**(-e), e
end subroutine

integer pure function getexp(x) result(e)
real(dp), intent(in) :: x
if (abs(x) < tiny(1._dp)) then
    e = 0
    return
endif
e = int(log10(abs(x)))
if (e < 0) e = e - 1
end function

end program
