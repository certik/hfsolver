module optimize

! Optimization algorithms

use types, only: dp
use utils, only: stop_error
implicit none
private
public bisect, brent, bracket, parabola_vertex, linregress

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
    if (abs(fc) < tiny(1._dp)) return   ! We need to make sure f(c) is not zero below
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
real(dp), parameter :: mintol = 1e-3_dp*epsilon(xa), &
    golden_ratio = (1+sqrt(5._dp))/2, cg = 2 - golden_ratio
real(dp) :: tol1, tol2, tmp1, tmp2, deltax, dx_temp, rat, xmid, &
        a, b, p, u, v, w, x, fu, fv, fw, fx
integer :: iter
logical :: verbose_
verbose_ = .true.
if (present(verbose)) verbose_ = verbose
verbose_ = .true.

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
rat = 0
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
            rat = p / tmp2            ! if parabolic step is useful.
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
        if ((fu <= fw) .or. abs(w - x) < tiny(1._dp)) then
            v = w; w = u; fv = fw; fw = fu
        else if ((fu <= fv) .or. abs(v - x) < tiny(1._dp) .or. abs(v - w) < tiny(1._dp)) then
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

subroutine bracket(f, xa, xb, xc, fa, fb, fc, grow_limit, maxiter, verbose)
procedure(func) :: f
real(dp), intent(inout) :: xa, xb
real(dp), intent(out) :: xc, fa, fb, fc
real(dp), intent(in) :: grow_limit
integer, intent(in) :: maxiter
logical, intent(in), optional :: verbose
real(dp), parameter :: golden_ratio = (1+sqrt(5._dp))/2
real(dp) :: denom, dum, fw, tmp1, tmp2, val, w, wlim
integer :: iter
logical :: verbose_
verbose_ = .false.
if (present(verbose)) verbose_ = verbose

fa = f(xa)
fb = f(xb)
if (fa < fb) then                      ! Switch so fa > fb
    dum = xa; xa = xb; xb = dum
    dum = fa; fa = fb; fb = dum
end if
xc = xb + golden_ratio*(xb - xa)
fc = f(xc)
iter = 0
do while (fc < fb)
    tmp1 = (xb - xa)*(fb - fc)
    tmp2 = (xb - xc)*(fb - fa)
    val = tmp2 - tmp1
    if (abs(val) < tiny(1.0_dp)) then
        denom = 2*tiny(1.0_dp)
    else
        denom = 2*val
    end if
    w = xb - ((xb - xc)*tmp2 - (xb - xa)*tmp1) / denom
    wlim = xb + grow_limit*(xc - xb)
    if (iter > maxiter) call stop_error("Too many iterations.")
    iter = iter + 1
    if ((w - xc)*(xb - w) > 0) then
        fw = f(w)
        if (fw < fc) then
            xa = xb; xb = w; fa = fb; fb = fw
            return
        else if (fw > fb) then
            xc = w; fc = fw
            return
        end if
        w = xc + golden_ratio*(xc - xb)
        fw = f(w)
    else if ((w - wlim)*(wlim - xc) >= 0) then
        w = wlim
        fw = f(w)
    else if ((w - wlim)*(xc - w) > 0) then
        fw = f(w)
        if (fw < fc) then
            xb = xc; xc = w; w = xc + golden_ratio*(xc - xb)
            fb = fc; fc = fw; fw = f(w)
        end if
    else
        w = xc + golden_ratio*(xc - xb)
        fw = f(w)
    end if
    xa = xb; xb = xc; xc = w
    fa = fb; fb = fc; fc = fw
    if (verbose_) then
        print "(i2, ':  xa = ', f23.12, ' xb = ', f23.12, ' xc = ', f23.12)", &
            iter, xa, xb, xc
    end if
end do
end subroutine

subroutine parabola_vertex(x1, y1, x2, y2, x3, y3, xv, yv)
! Calculates the position (xv, yv) of the parabola vertex.
! The parabola is defined by 3 points (x1, y1), (x2, y2), (x3, y3).
real(dp), intent(in) :: x1, y1, x2, y2, x3, y3
real(dp), intent(out) :: xv, yv
real(dp) :: denom, A, B, C
denom = (x1 - x2) * (x1 - x3) * (x2 - x3)
A     = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom
B     = (x3**2 * (y1 - y2) + x2**2 * (y3 - y1) + x1**2 * (y2 - y3)) / denom
C     = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + &
            x1 * x2 * (x1 - x2) * y3) / denom
xv = -B / (2*A)
yv = C - B**2 / (4*A)
end subroutine

subroutine linregress(x, y, slope, intercept, r, stderr_slope, stderr_intercept)
! Calculates simple linear regression of (x, y)
real(dp), intent(in) :: x(:), y(:)
real(dp), intent(out) :: slope ! y = intercept + slope*x
real(dp), intent(out) :: intercept ! y = intercept + slope*x
real(dp), intent(out) :: r ! correlation coefficient (the coefficient of
                           ! determination is r^2)
real(dp), intent(out) :: stderr_slope ! standard error of the slope
real(dp), intent(out) :: stderr_intercept ! standard error of the intercept
real(dp) :: xmean, ymean, varx, covxy, vary, r_den, mse
integer :: N
N = size(x)
xmean = sum(x)/N
ymean = sum(y)/N

varx  = dot_product(x, x) - N*xmean**2    ! = dot_product(x-xmean, x-xmean)
covxy = dot_product(x, y) - N*xmean*ymean ! = dot_product(x-xmean, y-ymean)
vary  = dot_product(y, y) - N*ymean**2    ! = dot_product(y-ymean, y-ymean)

slope = covxy / varx
intercept = ymean - slope*xmean

r_den = sqrt(varx * vary)
if (abs(r_den) < tiny(1._dp)) then
    r = 0
else
    r = covxy / r_den
    ! Normalize to [-1, 1] in case of numerical error propagation
    if (r > 1) then
        r = 1
    else if (r < -1) then
        r = -1
    end if
endif
! 'mse' is a mean square error (the sum of squared residuals divided by number
! of model parameters), which can be calculated directly as:
!   mse = sum((y-(slope*x+intercept))**2) / (N-2)
! But equivalently it can also be calculated in a faster way as:
mse = (1-r**2) * vary / (N-2)
stderr_slope = sqrt(mse / varx)
stderr_intercept = sqrt(mse * (1._dp/N + xmean**2/varx))
end subroutine

end module
