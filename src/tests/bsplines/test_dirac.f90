program test_dirac
use types, only: dp
use bsplines, only: bspline, bspline_der, bspline_der2
use mesh, only: meshexp
use quadrature, only: gauss_pts, gauss_wts
use linalg, only: eigh
use utils, only: stop_error
implicit none

integer, parameter :: n = 30, k = 6, Nq=7
integer, parameter :: N_intervals = n-k+1
integer, parameter :: Nq_total = Nq*N_intervals, Nb=n-1 + 1
real(dp) :: t(n+k), rmin, rmax, a
real(dp) :: xiq(Nq), wtq(Nq), xa, xb, jac, x(Nq)
real(dp), allocatable :: xq(:), wq(:), hq(:)
real(dp), allocatable :: B(:,:), Bp(:,:), Bpp(:,:), Bi(:), Bj(:), Bip(:), Bjp(:)
real(dp), allocatable :: Bipp(:), Bjpp(:)
real(dp), allocatable :: Vq(:), Vqp(:)
real(dp), allocatable :: Am(:,:), Bm(:,:), sol(:,:), lam(:)
real(dp), allocatable :: solsP(:,:), solsQ(:,:)
real(dp) :: En, c, beta
integer :: i, j, kappa, Z, l, relat, u

allocate(Am(2*Nb,2*Nb), Bm(2*Nb,2*Nb), sol(2*Nb,2*Nb), lam(2*Nb))
allocate(xq(Nq_total), wq(Nq_total), hq(Nq_total))
allocate(B(Nq_total, Nb), Bp(Nq_total, Nb), Bpp(Nq_total, Nb))
allocate(Bi(Nq_total), Bj(Nq_total), Bip(Nq_total), Bjp(Nq_total), Vq(Nq_total))
allocate(Bipp(Nq_total), Bjpp(Nq_total), Vqp(Nq_total))
allocate(solsP(Nq_total,Nb))
allocate(solsQ(Nq_total,Nb))

rmin = 1e-15_dp
rmax = 10
rmin = 0
rmax = 10
a = 6e5
a = 1e4
Z = 83
kappa = 2
c = 137.03599907_dp
beta = sqrt(kappa**2-(Z/c)**2)
print *, beta

t(:k-1) = rmin
t(k:n+1) = meshexp(rmin, rmax, a, N_intervals)
!t(k:n+1) = -1
t(n+2:) = rmax
!do i = 0, N_intervals
!    t(i+k) = (rmax-rmin)*(real(i, dp)/n)**6+rmin
!end do

print *, "Constructing quadrature rule"
! Loop over knot spans (intervals), and constract a global quadrature rule
! Integrals of a function hq evaluated at the points xq are calculated using:
! sum(wq*hq)
xiq = gauss_pts(Nq)
wtq = gauss_wts(Nq)
do i = 1, n-k+1
    xa = t(i+k-1)
    xb = t(i+k)
    jac = (xb-xa)/2
    x = (xiq(:)+1) * jac + xa
    xq((i-1)*Nq+1:i*Nq) = x
    wq((i-1)*Nq+1:i*Nq) = wtq*jac
end do

print *, "Evaluating basis functions"
! Evaluate basis functions and their derivatives on quadrature grid
! Skip the first and last B-spline (that's why Nb=n-2).
do i = 1, Nb-1
    B(:,i)   = bspline     (t, i, k, xq)
    Bp(:,i)  = bspline_der (t, i, k, xq)
    Bpp(:,i) = bspline_der2(t, i, k, xq)
end do
B(:,Nb)   = xq**beta
Bp(:,Nb)  = beta*xq**(beta-1)
Bpp(:,Nb) = beta*(beta-1)*xq**(beta-2)

print *, "Assembly"
! Construct matrices A and B
do i = 1, Nb
    do j = 1, Nb
        Bi = B(:,i)
        Bj = B(:,j)
        Bip = Bp(:,i)
        Bjp = Bp(:,j)
        Bipp = Bpp(:,i)
        Bjpp = Bpp(:,j)
        Vq = -Z/xq
        Vqp = Z/xq**2
        ! A11
        hq = -c**2*(Bi*xq**2*Bjpp +2*(Bi*xq*Bjp)-kappa*(kappa+1)*Bi*Bj) &
            +c**4*Bi*xq**2*Bj+Bi*xq**2*Vq**2*Bj+2*c**2*Bi*xq**2*Vq*Bj
        Am(i,j) = sum(wq*hq)
        ! A12
        hq = -c*(2*Bi*xq**2*Vq*Bjp+2*(1-kappa)*Bi*xq*Vq*Bj+Bi*xq**2*Vqp*Bj)
        Am(i,j+Nb) = sum(wq*hq)
        ! A21
        hq = c*(2*Bi*xq**2*Vq*Bjp+2*(1+kappa)*Bi*xq*Vq*Bj+Bi*xq**2*Vqp*Bj)
        Am(i+Nb,j) = sum(wq*hq)
        ! A22
        hq = -c**2*(Bi*xq**2*Bjpp +2*(Bi*xq*Bjp)-(-kappa)*(-kappa+1)*Bi*Bj) &
            +c**4*Bi*xq**2*Bj+Bi*xq**2*Vq**2*Bj-2*c**2*Bi*xq**2*Vq*Bj
        Am(i+Nb,j+Nb) = sum(wq*hq)

        ! B11
        hq = B(:,i)*B(:,j)*xq**2
        Bm(i,j) = sum(wq*hq)
        ! B22
        Bm(i+Nb,j+Nb) = Bm(i,j)
    end do
end do

print *, "Checking symmetry"
do j = 1, 2*Nb
    do i = 1, j-1
        if (max(abs(Am(i,j)), abs(Am(j,i))) > tiny(1._dp)) then
            if (abs(Am(i,j)-Am(j,i)) / max(abs(Am(i,j)), abs(Am(j,i))) &
                    > 1e-8_dp) then
                print *, i, j, Am(i,j)-Am(j,i), Am(i,j), Am(j,i)
                call stop_error("Am not symmetric")
            end if
        end if
        if (abs(Bm(i,j)-Bm(j,i)) > 1e-12_dp) call stop_error("Bm not symmetric")
   end do
end do


print *, "Eigensolver"
! Solve an eigenproblem
call eigh(Am, Bm, lam, sol)

print *, "n, energy, exact energy, error"
do i = 1, Nb
    if (kappa > 0) then
        l = kappa
        relat = 3
    else
        l = -kappa-1
        relat = 2
    end if
    En = E_nl(c, l+i, l, real(Z, dp), relat)
    lam(i) = sqrt(lam(i)) - c**2
    print "(i4, f30.8, f18.8, es12.2)", i, lam(i), En, abs(lam(i)-En)
end do

do i = 1, Nb
    solsP(:,i) = 0
    solsQ(:,i) = 0
    do j = 1, Nb
        solsP(:,i) = solsP(:,i) + sol(j,i)*B(:,j)
        solsQ(:,i) = solsQ(:,i) + sol(j+Nb,i)*B(:,j)
    end do
end do

open(newunit=u, file="dirac_sol.txt", status="replace")
write(u,*) xq
do i = 1, Nb
    write(u,*) solsP(:,i)
end do
close(u)

contains

real(dp) function E_nl(c, n, l, Z, relat)
! Calculates exact energy for the radial Schroedinger/Dirac equations
real(dp), intent(in) :: c, Z ! speed of light in atomic units
integer, intent(in) :: n, l, relat
! quantum numbers (n, l), atomic number (z)
! relat == 0 ... Schroedinger equation
! relat == 2 ... Dirac equation, spin up
! relat == 3 ... Dirac equation, spin down

integer :: skappa
real(dp) :: beta
if (.not. (l >= 0)) call stop_error("'l' must be positive or zero")
if (.not. (n > l)) call stop_error("'n' must be greater than 'l'")
if (l == 0 .and. relat == 3) call stop_error("Spin must be up for l==0.")
if (relat == 0) then
    E_nl = - Z**2 / (2.0_dp * n**2)
else
    if (relat == 2) then
        skappa = -l - 1
    else
        skappa = -l
    end if
    beta = sqrt(skappa**2 - Z**2 / c**2)
    E_nl = c**2/sqrt(1+Z**2/(n + skappa + beta)**2/c**2) - c**2
end if
end function

end program
