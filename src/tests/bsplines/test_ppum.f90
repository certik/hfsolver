program test_ppum
use types, only: dp
use bsplines, only: bspline, bspline_der, bspline_der2
use mesh, only: meshexp
use quadrature, only: gauss_pts, gauss_wts
use linalg, only: eigh
use utils, only: stop_error, str
implicit none

call do_ppum_basis(3, 0._dp, 1._dp, 3, 1.5_dp)

contains

    subroutine do_ppum_basis(p, xmin, xmax, Ne, alpha)
    integer, intent(in) :: p, Ne
    real(dp), intent(in) :: xmin, xmax, alpha
    integer :: Nq, N_intervals, Nq_total
    real(dp) :: rmin, rmax, dx
    real(dp) :: xa, xb, jac
    real(dp), allocatable :: xq(:), wq(:), hq(:), t(:), xiq(:), wtq(:), x(:)
    real(dp), allocatable :: B(:), Bp(:), Bpp(:)
    real(dp), allocatable :: mesh(:)
    integer :: i, u, bindex
    integer :: n, k
    k = p+1
    if (mod(p, 2) == 0) then
        n = k+2
    else
        n = k+1
    end if
    bindex = n/2+1
    print *, "Constructing B-spline basis n =", n, ", k =", k
    N_intervals = n-k+1

    allocate(mesh(2*Ne))
    dx = (xmax-xmin)/Ne*(alpha-1)/2
    mesh(1) = 0
    do i = 1, Ne-1
        mesh(2*i)   = (xmax-xmin)/Ne * i + xmin - dx
        mesh(2*i+1) = (xmax-xmin)/Ne * i + xmin + dx
    end do
    mesh(2*Ne) = xmax

    Nq = 64
    Nq_total = Nq*(2*Ne-1)

    allocate(t(n+k), xiq(Nq), wtq(Nq), x(Nq))
    allocate(xq(Nq_total), wq(Nq_total), hq(Nq_total))
    allocate(B(Nq_total), Bp(Nq_total), Bpp(Nq_total))

    ! Loop over the mesh, and constract a global quadrature rule. Integrals of a
    ! function hq evaluated at the points xq are calculated using: sum(wq*hq)
    xiq = gauss_pts(Nq)
    wtq = gauss_wts(Nq)
    do i = 1, 2*Ne-1
        xa = mesh(i)
        xb = mesh(i+1)
        jac = (xb-xa)/2
        x = (xiq(:)+1) * jac + xa
        xq((i-1)*Nq+1:i*Nq) = x
        wq((i-1)*Nq+1:i*Nq) = wtq*jac
    end do


    open(newunit=u, file="ppum.txt", status="replace")
    write(u,*) xq
    do i = 1, Ne
        dx = (xmax-xmin)/Ne*(alpha-1)/2
        rmin = (xmax-xmin)/Ne * (i-1) + xmin - dx
        rmax = (xmax-xmin)/Ne * i + xmin + dx
        t(:k-1) = rmin
        t(k:n+1) = meshexp(rmin, rmax, 1._dp, N_intervals)
        t(n+2:) = rmax
        B   = bspline     (t, bindex, k, xq)
        Bp  = bspline_der (t, bindex, k, xq)
        Bpp = bspline_der2(t, bindex, k, xq)
        write(u,*) B
    end do

    close(u)
    end subroutine

end program
