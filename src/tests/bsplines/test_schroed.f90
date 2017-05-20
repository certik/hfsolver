program test_schroed
use types, only: dp
use bsplines, only: bspline, bspline_der, bspline_der2
use mesh, only: meshexp
use quadrature, only: gauss_pts, gauss_wts
use linalg, only: eigh
use utils, only: stop_error
use schroed_util, only: lho, radial
implicit none

integer, parameter :: n = 30, k = 6, Nq=7
integer, parameter :: N_intervals = n-k+1
integer, parameter :: Nq_total = Nq*N_intervals, Nb=n-2
real(dp) :: t(n+k), rmin, rmax, a
real(dp) :: xiq(Nq), wtq(Nq), xa, xb, jac, x(Nq)
real(dp), allocatable :: xq(:), wq(:), hq(:)
real(dp), allocatable :: B(:,:), Bp(:,:), Bpp(:,:)
real(dp), allocatable :: Am(:,:), Bm(:,:), c(:,:), lam(:)
real(dp) :: En, condA, condB
integer :: i, j, u, l, Z

allocate(Am(Nb,Nb), Bm(Nb,Nb), c(Nb,Nb), lam(Nb))
allocate(xq(Nq_total), wq(Nq_total), hq(Nq_total))
allocate(B(Nq_total, n), Bp(Nq_total, n), Bpp(Nq_total, n))

rmin = -10
rmax = 10
a = 1
call construct_basis()
call lho(xq, wq, B(:,:Nb), Bp(:,:Nb), lam, condA, condB)


rmin = 0
rmax = 30
a = 30
Z = 2
l = 2
call construct_basis()
call radial(Nb, xq, wq, B, Bp, Z, l)

open(newunit=u, file="bspline_basis.txt", status="replace")
write(u,*) t
write(u,*) xq
do i = 1, n
    write(u,*) B(:,i)
end do
close(u)

contains

    subroutine construct_basis()
    t(:k-1) = rmin
    t(k:n+1) = meshexp(rmin, rmax, a, N_intervals)
    t(n+2:) = rmax

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
        B(:,i)   = bspline     (t, i+1, k, xq)
        Bp(:,i)  = bspline_der (t, i+1, k, xq)
        Bpp(:,i) = bspline_der2(t, i+1, k, xq)
    end do
    end subroutine

end program
