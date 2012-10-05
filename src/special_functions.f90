module special_functions
! Special functions
use types, only: dp
use constants, only: i_, pi
use utils, only: stop_error, str
implicit none
private
public wigner3j, getgaunt, getgauntr, Fm, Inu_asympt_sum, Inu_series, &
    Knu_asympt_sum, Inu_formula, Knu_formula, Inu_formula2, Knu_formula2

contains

real(dp) function fact(n) result(r)
integer, intent(in) :: n
r = gamma(n+1.0_dp)
end function

real(dp) function wigner3j(j1, j2, j3, m1, m2, m3) result(r)
! Calculates wigner3j symbols.
integer, intent(in) :: j1, j2, j3 ! angular momentum quantum numbers
integer, intent(in) :: m1, m2, m3 ! magnetic quantum numbers
integer :: k, k1, k2, l1, l2, l3, n1, n2, sgn
if ((abs(m1) > j1) .or. (abs(m2) > j2) .or. (abs(m3) > j3)) then
    print *, "Arguments:"
    print *, j1, j2, j3
    print *, m1, m2, m3
    call stop_error("wigner3j: invalid arguments")
end if
if ((j1 == 0) .and. (j2 == 0) .and. (j3 == 0)) then
    r = 1
    return
end if

l1 = j2 - j1 + j3
l2 = j1 - j2 + j3
l3 = j1 + j2 - j3
if ((m1 + m2 + m3 /= 0) .or. (l1 < 0) .or. (l2 < 0) .or. (l3 < 0)) then
    r = 0
    return
end if
n1 = j1 - m1
n2 = j2 + m2
k1 = max(0, j2-j3-m1, j1-j3+m2)
k2 = min(l3, n1, n2)
if (mod(k1 + j1 + j2 + m3, 2) /= 0) then
    sgn = -1
else
    sgn = 1
end if
r = 0
do k = k1, k2
    r = r + sgn * fr(l1, l1-n2+k) * fr(l2, l2-n1+k) * fr(l3, l3-k) / &
            (fact(k) * fact(n1-k) * fact(n2-k))
    sgn = -sgn
end do
r = r * sqrt(fr(j1+m1, l1) * fr(j2+m2, l2) * fr(j3+m3, l3) * &
        fr(j3-m3, 1+j1+j2+j3) * fact(j1-m1) * fact(j2-m2))
end function

real(dp) function fr(a, b) result(r)
! Returns a!/b!, for a,b >= 0
integer, intent(in) :: a, b
integer :: i
if (a < 0 .or. b < 0) call stop_error("fr(a,b): must satisfy a, b >= 0")
r = 1
if (a < b) then
    do i = a + 1, b
        r = r / i
    end do
else if (a > b) then
    do i = b + 1, a
        r = r * i
    end do
end if
end function

subroutine getgaunt(Lmax, ck)
! Calculates the Gaunt c^k(l1, m1, l2, m2) coefficients
! It returns an array of coefficients for all combinations of k, l1, m1, l2, m2:
!   ck(k, l1, m1, l2, m2) = c^k(l1, m1, l2, m2)
! Indices out of bounds mean that the coefficient is zero.
integer, intent(in) :: Lmax  ! max(l1, l2)
real(dp), allocatable, intent(out) :: ck(:, :, :, :, :) ! ck(k, l1, m1, l2, m2)
integer :: k, l1, m1, l2, m2
allocate(ck(0:2*Lmax, 0:Lmax, -Lmax:Lmax, 0:Lmax, -Lmax:Lmax))
ck = 0
do l1 = 0, Lmax
    do l2 = 0, Lmax
        do m1 = -l1, l1
            do m2 = -l2, l2
                do k = abs(l1-l2), l1+l2, 2
                    if (abs(m1-m2) > k) cycle
                    ck(k, l1, m1, l2, m2) = (-1)**(-m1) * &
                        sqrt(1._dp*(2*l1+1)*(2*l2+1)) * &
                        wigner3j(l1, k, l2, 0, 0, 0) * &
                        wigner3j(l1, k, l2, -m1, m1-m2, m2)
                end do
            end do
        end do
    end do
end do
end subroutine

subroutine getgauntr(Lmax, gr)
integer, intent(in) :: Lmax  ! max(l1, l2)
! gr(i, j, k, l, m, n) = <ij|kl|mn>_R
real(dp), allocatable, intent(out) :: gr(:, :, :, :, :, :)
real(dp), allocatable :: ck(:, :, :, :, :)
integer :: l1, m1, l2, m2, l3, m3
integer :: Lmax2
Lmax2 = 4*Lmax
allocate(ck(0:2*Lmax2, 0:Lmax2, -Lmax2:Lmax2, 0:Lmax2, -Lmax2:Lmax2))
call getgaunt(Lmax2, ck)
allocate(gr(0:Lmax, -Lmax:Lmax, 0:2*Lmax, -2*Lmax:2*Lmax, 0:Lmax, -Lmax:Lmax))
gr = 0
do l1 = 0, Lmax
    do l2 = 0, 2*Lmax
        do l3 = 0, Lmax
            do m1 = -l1, l1
                do m2 = -l2, l2
                    do m3 = -l3, l3
                        gr(l1, m1, l2, m2, l3, m3) = &
                            gauntr(l1, m1, l2, m2, l3, m3, ck, Lmax2)
                    end do
                end do
            end do
        end do
    end do
end do
end subroutine

real(dp) recursive function gauntr(l1, m1, l2, m2, l3, m3, ck, Lmax) result(r)
integer, intent(in) :: l1, m1, l2, m2, l3, m3, Lmax
real(dp), intent(in) :: ck(0:, 0:, -Lmax:, 0:, -Lmax:)
if (m1 /= 0 .and. m2 /= 0 .and. m3 /= 0) then
    ! Case A
    r = 2*ck(l2, l1, m2+m3, l3, m3) &
            * real(U(m1, m2+m3) * U(m2, m2) * U(m3, m3), dp) &
        + 2*ck(l2, l1, m2-m3, l3, -m3) &
            * real(U(m1, m2-m3) * U(m2, m2) * U(m3, -m3), dp)
else if (m1 /= 0 .and. m2 /= 0 .and. m3 == 0) then
    ! Case B
    r = 2*ck(l2, l1, m2, l3, 0) * real(U(m1, m2) * U(m2, m2), dp)
else if (m1 /= 0 .and. m2 == 0 .and. m3 /= 0) then
    ! Convert to case B
    r = gauntr(l1, m1, l3, m3, l2, m2, ck, Lmax)
else if (m1 == 0 .and. m2 /= 0 .and. m3 /= 0) then
    ! Convert to case B
    r = gauntr(l3, m3, l2, m2, l1, m1, ck, Lmax)
else if (m2 == 0 .and. m3 == 0) then
    ! Case C
    r = ck(l2, l1, 0, l3, 0) * delta(m1, 0)
else if (m1 == 0 .and. m3 == 0) then
    ! Convert to case C
    r = gauntr(l2, m2, l1, m1, l3, m3, ck, Lmax)
else if (m1 == 0 .and. m2 == 0) then
    ! Convert to case C
    r = gauntr(l3, m3, l1, m1, l2, m2, ck, Lmax)
else
    call stop_error("Internal error.")
end if

contains

integer function delta(a, b) result(r)
integer, intent(in) :: a, b
if (a == b) then
    r = 1
else
    r = 0
end if
end function

complex(dp) function U(mu, m) result(r)
integer, intent(in) :: mu, m
if (mu > 0) then
    r = (delta(mu, m) + (-1)**m * delta(mu, -m))/sqrt(2._dp)
else if (mu == 0) then
    r = delta(0, m)
else
    r = (delta(mu, -m) - (-1)**m * delta(mu, m))/(i_ * sqrt(2._dp))
end if
end function

end function


subroutine Fm(maxm, t, F)
! Calculates F_m(t) for m=0,1,..,maxm, where
!
!     F_m(t) = \int_0^1 u^(2m) e^(-tu^2) du
!
! And assigns the result to the array F(m) = F_m(t).
!
! Conditions on max, t: 0 <= maxm <= 500, t >= 0.
!
! This routine is tested for absolute accuracy 1e-15 and relative accuracy
! 1e-12 for all maxm = 0..500 and all values "t" from the interval
! 0 <= t <= 2e8.
!
! The algorithm is based on [1], all equations are references from there. The
! idea is to use series expansion for F_maxm(t) and then the recursive relation
! (24) downwards to calculate F_m(t) for m < maxm. For t >= maxm + 0.5, the
! series would take too many iterations to converge and also the downwards
! relation (24) becomes inaccurate, so we calculate F_0(t) directly and use
! (24) upwards.
!
! [1] I. Shavitt: Methods in Computational Physics (Academic Press Inc., New
! York, 1963), vol. 2
integer, intent(in) :: maxm
real(dp), intent(in) :: t
! The array F(0:maxm) will be equal to F(m) = F_m(t) for m = 0..maxm
real(dp), intent(out) :: F(0:)

real(dp) :: s, term
integer :: m
if (maxm < 0 .or. maxm > 500) &
    call stop_error("Fm: only works for 0 <= m <= 500")
if (t < 0) call stop_error("Fm: only works for t >= 0")
if (ubound(F, 1) /= maxm) call stop_error("Fm: invalid bounds on F")

if (t < maxm + 0.5_dp) then
    ! Series expansion for F_m(t), between equations (24) and (25).
    ! Since (2*t)/(2*maxm+1) < 1, this will converge fast:
    term = 1._dp / (2*maxm + 1)
    s = term
    m = 1
    do while (term/s > epsilon(1._dp))
        term = term * (2*t) / (2*maxm + 2 * m + 1)
        ! "s" will only change if term/s > machine eps
        s = s + term
        m = m + 1
    end do
    F(maxm) = s * exp(-t)
    ! Eq. (24) downwards, for t < maxm+0.5, this converges well:
    do m = maxm-1, 0, -1
        F(m) = (2*t*F(m + 1) + exp(-t)) / (2*m + 1)
    end do
else
    ! Eq. for F_0(t) on page 7:
    F(0) = 0.5_dp * sqrt(pi/t) * erf(sqrt(t))
    ! Eq. (24) upwards, for t >= maxm+0.5, this converges well:
    do m = 0, maxm-1
        F(m + 1) = ((2*m + 1)*F(m) - exp(-t)) / (2*t)
    end do
endif
end subroutine


real(dp) function Knu_asympt_sum(nu, x) result(s)
! If nu = k + 1/2 for k = 0, 1, 2, ..., then the series terminates (the result
! is a polynomial) and this function returns an exact result for all "x".
real(dp), intent(in) :: nu, x
real(dp) :: term, max_term
integer :: k
term = 1
s = term
max_term = abs(term)
k = 1
do while (term/s > epsilon(1._dp))
    term = term * (4*nu**2-(2*k-1)**2)/(8*k*x)
    if (abs(term) > max_term) max_term = abs(term)
    if (max_term > 1e100_dp) then
        call stop_error("Knu asymptotic series does not converge.")
    end if
    s = s + term
    k = k + 1
end do

if (max_term / abs(s) > 10) then
    call stop_error("Knu asymptotic series lost too many significant digits.")
end if
end function

real(dp) function Inu_asympt_sum(nu, x) result(s)
! If nu = k + 1/2 for k = 0, 1, 2, ..., then the series terminates (the result
! is a polynomial), however, unlike Knu_asympt_sum(), this polynomial is only
! valid asymptotically. The exact result contains functions sinh(x) and
! cosh(x). As such, sinh(x)/exp(x) = cosh(x)/exp(x) = 1/2 for
! x > -log(epsilon(1._dp))/2. Then this function returns exact result for the
! given floating point precision for all x satisfying the inequality.
real(dp), intent(in) :: nu, x
real(dp) :: term, max_term
integer :: k
term = 1
s = term
max_term = abs(term)
k = 1
do while (term/s > epsilon(1._dp))
    term = - term * (4*nu**2-(2*k-1)**2)/(8*k*x)
    if (abs(term) > max_term) max_term = abs(term)
    if (max_term > 1e100_dp) then
        call stop_error("Inu asymptotic series does not converge.")
    end if
    s = s + term
    k = k + 1
end do

if (max_term / abs(s) > 10) then
    call stop_error("Inu asymptotic series lost too many significant digits.")
end if
end function

real(dp) function hyp0f1(a, x) result(s)
real(dp), intent(in) :: a, x
real(dp) :: term, max_term
integer :: k
term = 1
s = term
max_term = abs(term)
k = 1
do while (term/s > epsilon(1._dp))
    term = term * x/(k*(a+k-1))
    if (abs(term) > max_term) max_term = abs(term)
    if (max_term > 1e100_dp) then
        call stop_error("hyp0f1: series does not converge.")
    end if
    s = s + term
    k = k + 1
end do

if (max_term / abs(s) > 10) then
    call stop_error("hyp0f1: series lost too many significant digits.")
end if
end function


real(dp) function Inu_series(nu, x) result(r)
! Returns the modified Bessel function of the first kind using the series
! expansion around x=0.
real(dp), intent(in) :: nu, x
r = (x/2)**nu / gamma(nu+1) * hyp0f1(nu+1, x**2/4)
end function

subroutine Knu_formula(maxk, x, K)
! Returns K(k) = K_{k+1/2}(x) / exp(-x) for k=-1,0,1,...,maxk.
integer, intent(in) :: maxk
real(dp), intent(in) :: x
real(dp), intent(out) :: K(-1:)
integer :: l
K(-1) = sqrt(pi/(2*x))
K(0) = K(-1)
do l = 0, maxk-1
    ! K_{nu+1} = K_{nu-1} + 2*nu*K_nu/x
    K(l + 1) = K(l-1) + (2*l + 1)*K(l)/x
end do
end subroutine

subroutine Inu_formula(maxk, x, I)
! Returns I(k) = I_{k+1/2}(x) / exp(x) for k=-1,0,1,...,maxk.
integer, intent(in) :: maxk
real(dp), intent(in) :: x
real(dp), intent(out) :: I(-1:)
integer :: l
I(-1) = (1+exp(-2*x)) / sqrt(2*pi*x)
! I(0) = (1-exp(-2*x)) / sqrt(2*pi*x)
if (x > -log(epsilon(1._dp))/2) then
    ! (1-exp(-2*x)) = 1 in floating point precission
    I(0) = 1 / sqrt(2*pi*x)
else
    ! (1-exp(-2*x)) = 2*sinh(x)/exp(x) to avoid cancellation
    I(0) = sqrt(2/(pi*x)) * sinh(x)/exp(x)
end if
do l = 0, maxk-1
    ! I_{nu+1} = I_{nu-1} - 2*nu*I_nu/x
    I(l + 1) = I(l-1) - (2*l + 1)*I(l)/x
end do
end subroutine

real(dp) elemental function esinh(x) result(r)
real(dp), intent(in) :: x
r = sinh(x) / exp(x)
end function

real(dp) elemental function ecosh(x) result(r)
real(dp), intent(in) :: x
r = cosh(x) / exp(x)
end function

real(dp) function Inu_formula2(k, x) result(r)
integer, intent(in) :: k
real(dp), intent(in) :: x
select case (k)
    case (0)
        r = esinh(x)
    case (1)
        if (x > 0.9_dp) then
            r = -esinh(x)/x + ecosh(x)
        else
            r = x**2/3 + x**4/30 + x**6/840 + x**8/45360 + x**10/3991680 + &
                x**12/518918400 + x**14/93405312e3_dp
            r = r / exp(x)
        end if
    case (2)
        if (x > 0.9_dp) then
            r = (3/x**2 + 1)*esinh(x) - 3*ecosh(x)/x
        else
            r = x**3/15 + x**5/210 + x**7/7560 + x**9/498960 + &
                x**11/51891840 + x**13/7783776e3_dp
            r = r / exp(x)
        end if
    case (3)
        if (x > 0.9_dp) then
            r = -(15/x**3 + 6/x)*esinh(x) + (15/x**2 + 1)*ecosh(x)
        else
            r = x**4/105 + x**6/1890 + x**8/83160 + x**10/6486480 + &
                x**12/778377600 + x**14/132324192e3_dp
            r = r / exp(x)
        end if
    case (4)
        if (x < 0.2_dp) then
            r = x**5/945 + x**7/20790 + x**9/1081080 + x**11/97297200 + &
                    x**13/132324192e2_dp
            r = r/exp(x)
        else if (x < 1.7_dp) then
            r =  (5.1962956007054264229527112912e-17_dp -  &
                    1.7091983797257865101834560981e-15_dp*x +  &
                    2.54701673734823215056998758576e-14_dp*x**2 -  &
                    2.28250221664547749050764003136e-13_dp*x**3 +  &
                    1.37741730213761458713328353583e-12_dp*x**4 +  &
                    0.00105820105225973559701867830214_dp*x**5 +  &
                    0.00061204753548652964602354522028_dp*x**6 +  &
                    0.0000122798136500667925710228809483_dp*x**7 +  &
                    0.0000103551941240484592046206146607_dp*x**8)/ &
                (1 + 0.578384903100157953572414461702_dp*x -  &
                    0.0338500781802291050601812074314_dp*x**2 -  &
                    0.0165046449338415288740437462715_dp*x**3 +  &
                    0.000664629591462297572628738277521_dp*x**4 +  &
                    0.000244501438485835213668147591168_dp*x**5 -  &
                    0.0000102233345389484523128091507614_dp*x**6 -  &
                    2.37680351471245713258958333627e-6_dp*x**7 +  &
                    1.7722924118248819944246867464e-7_dp*x**8)
            r = r/exp(x)
        else if (x < 4) then
            r = (-1.53382275599031353712671171402e-7_dp + &
                    1.02445790975414106564424700453e-6_dp*x - &
                    3.20332766010782394861350740562e-6_dp*x**2 + &
                    6.22984887800335387305228161855e-6_dp*x**3 - &
                    8.44422432663036167543192411535e-6_dp*x**4 + &
                    0.00106667316609330486086491285357_dp*x**5 + &
                    0.000296879521931548405158141021032_dp*x**6 + &
                    9.16131185419642061414089694738e-6_dp*x**7 + &
                    4.84398673121860851791923295866e-6_dp*x**8)/ &
                (1 + 0.286715277505499007952029608584_dp*x - &
                    0.0405239647887260390437917379491_dp*x**2 - &
                    0.00665840458518766453585847552188_dp*x**3 + &
                    0.000270261982447901375808140464959_dp*x**4 + &
                    0.000271851192449937947349585685247_dp*x**5 - &
                    0.0000431403377516821653299354827079_dp*x**6 + &
                    2.65420161107604023280035849057e-6_dp*x**7 - &
                    6.20361275935376828687251028769e-8_dp*x**8)
            r = r/exp(x)
        else if (x < 10) then
            ! Produced by:
            ! FortranForm[MiniMaxApproximation[((105/x^4 + 45/x^2 +
            !   1)*Sinh[x]-(105/x^3 + 10/x)*Cosh[x])/Exp[x], {x, {4, 10}, 8, 8},
            !   WorkingPrecision->30]]
            r = (0.000395502959013236968661582656143_dp - &
                    0.001434648369704841686633794071_dp*x + &
                    0.00248783474583503473135143644434_dp*x**2 - &
                    0.00274477921388295929464613063609_dp*x**3 + &
                    0.00216275018107657273725589740499_dp*x**4 - &
                    0.000236779926184242197820134964535_dp*x**5 + &
                    0.0000882030507076791807159699814428_dp*x**6 - &
                    4.62078105288798755556136693122e-6_dp*x**7 + &
                    8.23671374777791529292655504214e-7_dp*x**8)/ &
                (1 + 0.504839286873735708062045336271_dp*x + &
                    0.176683950009401712892997268723_dp*x**2 + &
                    0.0438594911840609324095487447279_dp*x**3 + &
                    0.00829753062428409331123592322788_dp*x**4 + &
                    0.00111693697900468156881720995034_dp*x**5 + &
                    0.000174719963536517752971223459247_dp*x**6 + &
                    7.22885338737473776714257581233e-6_dp*x**7 + &
                    1.64737453771748367647332279826e-6_dp*x**8)
        else if (x < 20) then
            r = (105/x**4 + 45/x**2 + 1)*esinh(x) - (105/x**3 + 10/x)*ecosh(x)
        else
            r = (105/x**4 - 105/x**3 + 45/x**2 - 10/x + 1)/2
        end if
    case default
        call stop_error("k = " // str(k) // " not implemented.")
end select
r = r * sqrt(2/(pi*x))
end function

real(dp) function Knu_formula2(k, x) result(r)
integer, intent(in) :: k
real(dp), intent(in) :: x
select case (k)
    case (0)
        r = 1
    case (1)
        r = 1/x + 1
    case (2)
        r = 3/x**2 + 3/x + 1
    case (3)
        r = 15/x**3 + 15/x**2 + 6/x + 1
    case (4)
        r = 105/x**4 + 105/x**3 + 45/x**2 + 10/x + 1
    case default
        call stop_error("k = " // str(k) // " not implemented.")
end select
r = r * sqrt(pi/(2*x))
end function

end module
