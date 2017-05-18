program test_schroed
use types, only: dp
use bsplines, only: bspline, bspline_der, bspline_der2
use mesh, only: meshexp
use quadrature, only: gauss_pts, gauss_wts
use linalg, only: eigh
use utils, only: stop_error
implicit none

integer, parameter :: n = 40, k = 6, Nq=64
integer, parameter :: N_intervals = n-k+1
integer, parameter :: Nq_total = Nq*N_intervals, Nb=n-2
real(dp) :: t(n+k), rmin, rmax, a
real(dp) :: xiq(Nq), wtq(Nq), xa, xb, jac, x(Nq)
real(dp), allocatable :: xq(:), wq(:), hq(:)
real(dp), allocatable :: B(:,:), Bp(:,:), Bpp(:,:)
real(dp), allocatable :: Am(:,:), Bm(:,:), c(:,:), lam(:)
real(dp) :: En
integer :: i, j, u, l, Z

allocate(Am(Nb,Nb), Bm(Nb,Nb), c(Nb,Nb), lam(Nb))
allocate(xq(Nq_total), wq(Nq_total), hq(Nq_total))
allocate(B(Nq_total, n), Bp(Nq_total, n), Bpp(Nq_total, n))

rmin = -10
rmax = 10
a = 1
Z = 2
l = 2

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

print *, "Assembly"
! Construct matrices A and B
do i = 1, Nb
    do j = 1, Nb
        ! A
        ! Both of these work:
        !hq = -B(:,i)*Bpp(:,j)/2 + B(:,i)*B(:,j)*(-Z/xq+l*(l+1)/(2*xq**2))
        hq = Bp(:,i)*Bp(:,j)/2 + B(:,i)*B(:,j)*(xq**2)/2
        Am(i,j) = sum(wq*hq)

        ! B
        hq = B(:,i)*B(:,j)
        Bm(i,j) = sum(wq*hq)
    end do
end do

print *, "Checking symmetry"
do j = 1, Nb
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
call eigh(Am, Bm, lam, c)

print *, "n, energy, exact energy, error"
do i = 1, Nb
    En = 0.5_dp + (i-1)
    print "(i4, f30.8, f18.8, es12.2)", i, lam(i), En, abs(lam(i)-En)
end do

open(newunit=u, file="bspline_basis.txt", status="replace")
write(u,*) t
write(u,*) xq
do i = 1, n
    write(u,*) B(:,i)
end do
close(u)

end program
