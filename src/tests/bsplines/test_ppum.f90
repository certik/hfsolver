program test_ppum
use types, only: dp
use bsplines, only: bspline, bspline_der, bspline_der2
use mesh, only: meshexp
use quadrature, only: gauss_pts, gauss_wts
use linalg, only: eigh
use utils, only: stop_error, str
implicit none

call do_ppum_basis(3)

contains

    subroutine do_ppum_basis(p)
    integer, intent(in) :: p
    integer :: Nq, N_intervals, Nq_total, Nb
    real(dp) :: rmin, rmax, a
    real(dp) :: xa, xb, jac
    real(dp), allocatable :: xq(:), wq(:), hq(:), t(:), xiq(:), wtq(:), x(:)
    real(dp), allocatable :: B(:,:), Bp(:,:), Bpp(:,:)
    integer :: i, u
    integer :: n, k
    k = p+1
    n = k+1
    print *, "Constructing B-spline basis n =", n, ", k =", k
    Nq = 64
    N_intervals = n-k+1
    Nq_total = Nq*N_intervals
    Nb=n

    allocate(t(n+k), xiq(Nq), wtq(Nq), x(Nq))
    allocate(xq(Nq_total), wq(Nq_total), hq(Nq_total))
    allocate(B(Nq_total, n), Bp(Nq_total, n), Bpp(Nq_total, n))

    rmin = 0
    rmax = 1
    a = 1
    t(:k-1) = rmin
    t(k:n+1) = meshexp(rmin, rmax, a, N_intervals)
    t(n+2:) = rmax

    ! Loop over knot spans (intervals), and constract a global quadrature rule
    ! Integrals of a function hq evaluated at the points xq are calculated
    ! using: sum(wq*hq)
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


    open(newunit=u, file="ppum.txt", status="replace")
    write(u,*) xq

    rmin = 0
    rmax = 5._dp/12
    t(:k-1) = rmin
    t(k:n+1) = meshexp(rmin, rmax, a, N_intervals)
    t(n+2:) = rmax
    do i = 1, Nb
        B(:,i)   = bspline     (t, i, k, xq)
        Bp(:,i)  = bspline_der (t, i, k, xq)
        Bpp(:,i) = bspline_der2(t, i, k, xq)
    end do
    write(u,*) B(:,Nb/2+1)

    rmin = 1._dp/4
    rmax = 3._dp/4
    t(:k-1) = rmin
    t(k:n+1) = meshexp(rmin, rmax, a, N_intervals)
    t(n+2:) = rmax
    do i = 1, Nb
        B(:,i)   = bspline     (t, i, k, xq)
        Bp(:,i)  = bspline_der (t, i, k, xq)
        Bpp(:,i) = bspline_der2(t, i, k, xq)
    end do
    write(u,*) B(:,Nb/2+1)

    rmin = 7._dp/12
    rmax = 1
    t(:k-1) = rmin
    t(k:n+1) = meshexp(rmin, rmax, a, N_intervals)
    t(n+2:) = rmax
    do i = 1, Nb
        B(:,i)   = bspline     (t, i, k, xq)
        Bp(:,i)  = bspline_der (t, i, k, xq)
        Bpp(:,i) = bspline_der2(t, i, k, xq)
    end do
    write(u,*) B(:,Nb/2+1)

    close(u)
    end subroutine

end program
