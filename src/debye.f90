module debye
use types, only: dp
use debye_potential, only: V0, V1, V2, V3, V4, qp
use debye_potential_series, only: S0, S1, S2, S3, S4
use utils, only: stop_error, str
use constants, only: pi
use special_functions, only: Inu_asympt_sum, Inu_series, Knu_asympt_sum, &
    Inu_formula, Knu_formula, Inu_formula2
implicit none

private
public Vk, Sk, Vk2, Vk3

contains

real(dp) function Vk(k, D, r1, r2) result(res)
! Uses direct formula. Depending on r1, r2, typically the result cannot be
! trusted for D >= 100. The accuracy gets worse with increasing "k".
integer, intent(in) :: k
real(dp), intent(in) :: D, r1, r2
real(qp) :: D_, r1_, r2_, res_
D_ = real(D, qp)
r1_ = real(r1, qp)
r2_ = real(r2, qp)
select case (k)
    case (0)
        res_ = V0(D_, r1_, r2_)
    case (1)
        res_ = V1(D_, r1_, r2_)
    case (2)
        res_ = V2(D_, r1_, r2_)
    case (3)
        res_ = V3(D_, r1_, r2_)
    case (4)
        res_ = V4(D_, r1_, r2_)
    case default
        call stop_error("k = " // str(k) // " not implemented.")
end select
res = real(res_, dp)
end function

real(dp) function Sk(k, D, r1, r2) result(res)
! Uses series expansion in "alpha".
! The convergence depends on two parameters:
!     alpha = max(r1, r2) / D
!     t = min(r1, r2) / max(r1, r2)
! The series satisfies the folling properties:
!     I  : From the definition, t <= 1
!     II : The series contains a few "t" terms for each power of alpha
!     III: For a given power of alpha, all "t" terms have equal or higher power
! Three cases can happen in terms of convergence:
!     a) If alpha <= 1, the series always converges (I and II)
!     b) If alpha >  1, but alpha * t <= 1, the series always converges (III)
!     c) If alpha >  1  and alpha * t >  1, the series can diverge and in this
!        case the results are usually not accurate. A warning is issed.
! In the a) and b) cases, the accuracy should be at least 8 significant
! digits. In the c) case, typically the accuracy is much lower, sometimes not
! even a single digit is correct.
! There doesn't seem to be significant loss of accuracy with increasing "k".
integer, intent(in) :: k
real(dp), intent(in) :: D, r1, r2
real(dp) :: alpha, t, D_, r1_, r2_, res_, rmin, rmax
D_ = real(D, dp)
r1_ = real(r1, dp)
r2_ = real(r2, dp)
! r1 >= r2
rmin = r2_
rmax = r1_
alpha = rmax / D_
t = rmin / rmax
! Don't issue the warning, because we also need to allow t > 1.
!if (alpha > 1 .and. alpha * t > 1) then
!    print *, "alpha =", alpha, "t =", t
!    print *, "WARNING: alpha > 1 and alpha*t > 1; The series for V_k(r1, r2)"
!    print *, "         might diverge. The results cannot be trusted."
!end if
select case (k)
    case (0)
        res_ = S0(alpha, t)
    case (1)
        res_ = S1(alpha, t)
    case (2)
        res_ = S2(alpha, t)
    case (3)
        res_ = S3(alpha, t)
    case (4)
        res_ = S4(alpha, t)
    case default
        call stop_error("k = " // str(k) // " not implemented.")
end select
res_ = res_ * exp(-alpha) / rmax
res = real(res_, dp)
end function

real(dp) function Vk2(k, D, r1, r2) result(r)
! Uses the modified Bessel functions
integer, intent(in) :: k
real(dp), intent(in) :: D, r1, r2
real(dp) :: rmin, rmax, a, b, C
rmin = r2
rmax = r1
a = Knu_asympt_sum(k+0.5_dp, rmax/D) * sqrt(pi/(2*rmax/D))
if (rmin/D > -log(epsilon(1._dp))/2) then
    b = Inu_asympt_sum(k+0.5_dp, rmin/D) / sqrt(2*pi*rmin/D)
    C = exp((rmin-rmax)/D)
else
    a = a * exp(-rmax/D)
    b = Inu_series(k+0.5_dp, rmin/D)
    C = 1
end if
r = (2*k+1)*C*a*b/sqrt(rmin*rmax)
end function

real(dp) function Vk3(k, D, r1, r2) result(r)
! Uses the modified Bessel functions
integer, intent(in) :: k
real(dp), intent(in) :: D, r1, r2
real(dp) :: rmin, rmax, a, b, C
real(dp), dimension(-1:k) :: I_, K_
rmin = r2
rmax = r1
call Inu_formula(k, rmin/D, I_)
call Knu_formula(k, rmax/D, K_)
a = K_(k)
b = I_(k)
C = exp((rmin-rmax)/D)
r = (2*k+1)*C*a*b/sqrt(rmin*rmax)
end function

end module
