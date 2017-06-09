program test_dirac
use types, only: dp
use bsplines, only: bspline, bspline_der, bspline_der2
use mesh, only: meshexp
use quadrature, only: gauss_pts, gauss_wts
use linalg, only: eigh
use utils, only: stop_error
use radial_util, only: radial_dirac
implicit none

integer, parameter :: n = 30, k = 6, Nq=7
integer, parameter :: N_intervals = n-k+1
integer, parameter :: Nq_total = Nq*N_intervals, Nb=n-1
real(dp) :: t(n+k), rmin, rmax, a
real(dp) :: xiq(Nq), wtq(Nq), xa, xb, jac, x(Nq)
real(dp), allocatable :: xq(:), wq(:)
real(dp), allocatable :: B(:,:), Bp(:,:), Bpp(:,:)
real(dp) :: c
integer :: i, kappa, Z

allocate(xq(Nq_total), wq(Nq_total))
allocate(B(Nq_total, Nb), Bp(Nq_total, Nb), Bpp(Nq_total, Nb))

rmin = 1e-15_dp
rmax = 10
rmin = 0
rmax = 10
a = 6e5
a = 1e4
Z = 83
kappa = 2
c = 137.03599907_dp

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
do i = 1, Nb
    B(:,i)   = bspline     (t, i, k, xq)
    Bp(:,i)  = bspline_der (t, i, k, xq)
    Bpp(:,i) = bspline_der2(t, i, k, xq)
end do

call radial_dirac(Nb, xq, wq, B, Bp, Z, kappa, c)

end program
