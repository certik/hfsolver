module optimize

! Optimization algorithms

use types, only: dp
use utils, only: stop_error
implicit none
private
public bisect, brent

interface
    real(dp) function func(x)
    import :: dp
    implicit none
    real(dp), intent(in) :: x
    end function
end interface

contains

real(dp) function bisect(f, a, b, tol) result(c)
! Solves f(x) = 0 on the interval [a, b] using the bisection method
procedure(func) :: f
real(dp), intent(in) :: a, b, tol
real(dp) :: a_, b_, fa, fb, fc
a_ = a; b_ = b
fa = f(a_)
fb = f(b_)
if (fa * fb >= 0) then
    call stop_error("bisect: f(a) and f(b) must have opposite signs")
end if
do while (b_ - a_ > tol)
    c = (a_ + b_) / 2
    fc = f(c)
    if (fc == 0) return   ! We need to make sure f(c) is not zero below
    if (fa * fc < 0) then
        b_ = c
        fb = fc
    else
        a_ = c
        fa = fc
    end if
end do
c = (a_ + b_)/2
end function

subroutine brent(f, xa, xb, xc, tol, maxiter, xmin, fxmin, verbose)
procedure(func) :: f
real(dp), intent(in) :: xa, xb, xc, tol
integer, intent(in) :: maxiter
real(dp), intent(out) :: xmin, fxmin
logical, intent(in), optional :: verbose
real(dp), parameter :: mintol = 1e-12_dp, golden_ratio = (1+sqrt(5._dp))/2, &
        cg = 2 - golden_ratio
real(dp) :: tol1, tol2, tmp1, tmp2, deltax, dx_temp, rat, xmid, &
        a, b, p, u, v, w, x, fu, fv, fw, fx
integer :: iter
logical :: verbose_
verbose_ = .false.
if (present(verbose)) verbose_ = verbose

x = xb; w = xb; v = xb
fx = f(x)
fw = fx; fv = fx
if (xa < xc) then
    a = xa
    b = xc
else
    a = xc
    b = xa
end if
deltax = 0
do iter = 1, maxiter
    if (verbose_) then
        print "(i2, ':  x = ', f23.12, '     tol = ', es10.2)", iter, x, &
            (b - a) / 2
    end if
    tol1 = tol*abs(x) + mintol
    tol2 = 2*tol1
    xmid = 0.5_dp*(a + b)
    if (abs(x - xmid) < (tol2 - 0.5_dp*(b - a))) then ! check for convergence
        xmin = x
        fxmin = fx
        return
    end if
    if (abs(deltax) <= tol1) then
        if (x >= xmid) then           ! do a golden section step
            deltax = a - x
        else
            deltax = b - x
        end if
        rat = cg*deltax
    else
        tmp1 = (x - w)*(fx - fv)      ! do a parabolic step
        tmp2 = (x - v)*(fx - fw)
        p = (x - v)*tmp2 - (x - w)*tmp1;
        tmp2 = 2*(tmp2 - tmp1)
        if (tmp2 > 0) then
            p = -p
        end if
        tmp2 = abs(tmp2)
        dx_temp = deltax
        deltax = rat
        if ((p > tmp2*(a - x)) .and. (p < tmp2*(b - x)) .and. &
                (abs(p) < abs(0.5_dp*tmp2*dx_temp))) then ! check parabolic fit
            rat = p*1.0 / tmp2        ! if parabolic step is useful.
            u = x + rat
            if ((u - a) < tol2 .or. (b - u) < tol2) then
                if (xmid - x >= 0) then
                    rat = tol1
                else
                    rat = -tol1
                end if
            end if
        else
            if (x >= xmid) then       ! if it's not do a golden section step
                deltax = a - x
            else
                deltax = b - x
            end if
            rat = cg*deltax
        end if
    end if

    if (abs(rat) < tol1) then         ! update by at least tol1
        if (rat >= 0) then
            u = x + tol1
        else
            u = x - tol1
        end if
    else
        u = x + rat
    end if
    fu = f(u) ! calculate new output value

    if (fu > fx) then                 ! if it's bigger than current
        if (u < x) then
            a = u
        else
            b = u
        end if
        if ((fu <= fw) .or. (w == x)) then
            v = w; w = u; fv = fw; fw = fu
        else if ((fu <= fv) .or. (v == x) .or. (v == w)) then
            v = u; fv = fu
        end if
    else
        if (u >= x) then
            a = x
        else
            b = x
        end if
        v = w; w = x; x = u
        fv = fw; fw = fx; fx = fu
    end if
end do
call stop_error("brent: The maximum number of iterations exceeded.")
end subroutine


end module
