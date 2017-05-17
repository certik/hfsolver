program test_ppum
use types, only: dp
use bsplines, only: bspline, bspline_der, bspline_der2
use mesh, only: meshexp
use quadrature, only: gauss_pts, gauss_wts
use linalg, only: eigh
use utils, only: stop_error, str
implicit none

call do_ppum_basis(3, 0._dp, 1._dp, 4, 6, 1.5_dp)

contains

    subroutine do_ppum_basis(p, xmin, xmax, Ne, Nenr, alpha)
    integer, intent(in) :: p, Ne, Nenr
    real(dp), intent(in) :: xmin, xmax, alpha
    integer :: Nq, N_intervals, Nq_total
    real(dp) :: rmin, rmax, dx
    real(dp) :: xa, xb, jac, x0
    real(dp), allocatable :: xq(:), wq(:), hq(:), t(:), xiq(:), wtq(:), x(:)
    real(dp), allocatable :: W(:,:), Wp(:,:), Wpp(:,:), S(:), Sp(:), Spp(:)
    real(dp), allocatable :: wi(:,:), wip(:,:), wipp(:,:)
    real(dp), allocatable :: enr(:,:,:), enrp(:,:,:)
    real(dp), allocatable :: mesh(:)
    integer :: i, j, u, bindex
    integer :: n, k
    k = p+1
    if (mod(p, 2) == 0) then
        n = k+2
    else
        if (p == 1) then
            n = k+1
        else
            n = k+3
        end if
    end if
    bindex = n/2+1
    print *, "Constructing B-spline basis n =", n, ", k =", k
    N_intervals = n-k+1

    allocate(mesh(2*Ne))
    dx = (xmax-xmin)/Ne*(alpha-1)/2
    mesh(1) = xmin
    do i = 1, Ne-1
        mesh(2*i)   = (xmax-xmin)/Ne * i + xmin - dx
        mesh(2*i+1) = (xmax-xmin)/Ne * i + xmin + dx
    end do
    mesh(2*Ne) = xmax

    Nq = 64
    Nq_total = Nq*(2*Ne-1)

    allocate(t(n+k), xiq(Nq), wtq(Nq), x(Nq))
    allocate(xq(Nq_total), wq(Nq_total), hq(Nq_total))
    allocate(W(Nq_total,Ne), Wp(Nq_total,Ne), Wpp(Nq_total,Ne))
    allocate(S(Nq_total), wi(Nq_total,Ne))
    allocate(Sp(Nq_total), wip(Nq_total,Ne))
    allocate(Spp(Nq_total), wipp(Nq_total,Ne))
    allocate(enr(Nq_total,Nenr,Ne))
    allocate(enrp(Nq_total,Nenr,Ne))

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

    ! Construct weight functions Wi(x)
    do i = 1, Ne
        dx = (xmax-xmin)/Ne*(alpha-1)/2
        rmin = (xmax-xmin)/Ne * (i-1) + xmin - dx
        rmax = (xmax-xmin)/Ne * i + xmin + dx
        t(:k-1) = rmin
        t(k:n+1) = meshexp(rmin, rmax, 1._dp, N_intervals)
        t(n+2:) = rmax
        W(:,i)   = bspline     (t, bindex, k, xq)
        Wp(:,i)  = bspline_der (t, bindex, k, xq)
        Wpp(:,i) = bspline_der2(t, bindex, k, xq)
    end do

    ! Construct PU functions using Shepard constructions
    S = 0
    Sp = 0
    Spp = 0
    do i = 1, Ne
        S = S + W(:,i)
        Sp = Sp + Wp(:,i)
        Spp = Spp + Wpp(:,i)
    end do
    do i = 1, Ne
        wi(:,i) = W(:,i)/S
        wip(:,i) = (Wp(:,i)*S - W(:,i)*Sp)/S**2
        wipp(:,i) = (Wpp(:,i)*S**2 - 2*Wp(:,i)*S*Sp - W(:,i)*S*Spp &
            + 2*W(:,i)*Sp**2)/S**3
    end do

    ! Construct enrichment functions enr(x)
    do i = 1, Ne
        dx = (xmax-xmin)/Ne*(alpha-1)/2
        rmin = (xmax-xmin)/Ne * (i-1) + xmin - dx
        rmax = (xmax-xmin)/Ne * i + xmin + dx
        x0 = (rmin+rmax)/2
        jac = (rmax-rmin)/2
        do j = 1, Nenr
            enr(:,j,i) = legendre_p((xq-rmin)/jac-1, j-1)
            enrp(:,j,i) = legendre_p_der((xq-rmin)/jac-1, j-1)
            where (xq < rmin .or. xq > rmax)
                enr(:,j,i) = 0
                enrp(:,j,i) = 0
            end where
        end do
    end do


    open(newunit=u, file="ppum.txt", status="replace")
    write(u,*) xq
    do i = 1, Ne
        write(u,*) W(:,i)
    end do
    do i = 1, Ne
        write(u,*) wi(:,i)
    end do
    do i = 1, Ne
        write(u,*) wip(:,i)
    end do
    do i = 1, Ne
        write(u,*) wipp(:,i)
    end do
    close(u)

    open(newunit=u, file="enr.txt", status="replace")
    do i = 1, Ne
        do j = 1, Nenr
            write(u,*) enr(:,j,i)
            write(u,*) enrp(:,j,i)
        end do
    end do
    close(u)
    end subroutine

    function legendre_p(x, n) result(r)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: n
    real(dp) :: r(size(x))
    select case (n)
        case (0)
            r = 1
        case (1)
            r = x
        case (2)
            r = (3*x**2-1)/2
        case (3)
            r = (5*x**3-3*x)/2
        case (4)
            r = (35*x**4-30*x**2+3)/8
        case (5)
            r = (63*x**5-70*x**3+15*x)/8
        case (6)
            r = (231*x**6-315*x**4+105*x**2-5)/16
        case (7)
            r = (429*x**7-693*x**5+315*x**3-35*x)/16
        case (8)
            r = (6435*x**8-12012*x**6+6930*x**4-1260*x**2+35)/128
        case (9)
            r = (12155*x**9-25740*x**7+18018*x**5-4620*x**3+315*x)/128
        case (10)
            r = (46189*x**10-109395*x**8+90090*x**6-30030*x**4+3465*x**2-63)/256
        case default
            call stop_error("n is too high")
    end select
    end function

    function legendre_p_der(x, n) result(r)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: n
    real(dp) :: r(size(x))
    select case (n)
        case (0)
            r = 0
        case (1)
            r = 1
        case (2)
            r = (6*x)/2
        case (3)
            r = (15*x**2-3)/2
        case (4)
            r = (140*x**3-60*x)/8
        case (5)
            r = (315*x**4-210*x**2+15)/8
        case (6)
            r = (1386*x**5-1260*x**3+210*x)/16
        case (7)
            r = (3003*x**6-3465*x**4+945*x**2-35)/16
        case (8)
            r = (51480*x**7-72072*x**5+27720*x**3-2520*x)/128
        case (9)
            r = (109395*x**8-180180*x**6+90090*x**4-13860*x**2+315)/128
        case (10)
            r = (461890*x**9-875160*x**7+540540*x**5-120120*x**3+6930*x)/256
        case default
            call stop_error("n is too high")
    end select
    end function

end program