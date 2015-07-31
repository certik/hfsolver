module xc

! Exchange and correlation potentials

use types, only: dp
use constants, only: pi
implicit none
private
public get_Vxc_vwn, xc_pz

contains

subroutine get_Vxc_vwn(R, rho, relat, c, exc, Vxc)
real(dp), intent(in) :: R(:) ! radial grid
real(dp), intent(in) :: rho(:) ! charge density
logical, intent(in) :: relat ! .true. return RLDA, otherwise LDA
real(dp), intent(in) :: c ! speed of light
real(dp), intent(out) :: Vxc(:), exc(:)

integer :: i
do i = 1, size(R)
    call xc_vwn(rho(i), relat, c, exc(i), Vxc(i))
end do
end subroutine

subroutine xc_vwn(n, relat, c_light, exc, Vxc)
! Calculates XC LDA density and potential from the charge density "n".
real(dp), intent(in) :: n ! charge density (scalar)
real(dp), intent(in) :: c_light ! speed of light
logical, intent(in) :: relat ! if .true. returns RLDA, otherwise LDA
real(dp), intent(out) :: exc ! XC density
real(dp), intent(out) :: Vxc ! XC potential

real(dp), parameter :: y0 = -0.10498_dp
real(dp), parameter :: b = 3.72744_dp
real(dp), parameter :: c = 12.9352_dp
real(dp), parameter :: A = 0.0621814_dp

real(dp) :: Q, rs, y, ec, ex, Vc, Vx, beta, mu, R, S

if (abs(n) < tiny(1._dp)) then
    exc = 0
    Vxc = 0
    return
end if

Q = sqrt(4*c - b**2)
rs = (3/(4*pi*n))**(1.0_dp/3)
y = sqrt(rs)
ec = A/2 * (log(y**2/get_Y(y, b, c)) + 2*b/Q * atan(Q/(2*y+b))  &
   - b*y0/get_Y(y0, b, c) * ( &
            log((y-y0)**2 / get_Y(y, b, c)) &
            + 2*(b+2*y0) / Q * atan(Q/(2*y+b)) &
          ) )
Vc = ec - A/6 * (c*(y-y0)-b*y0*y)/((y-y0)*get_Y(y, b, c))
ex = -3/(4*pi) * (3*pi**2*n)**(1.0_dp/3)
Vx = 4*ex/3

if (relat) then
    beta = -4 * pi * ex / (3 * c_light)
    mu = sqrt(1 + beta**2)
    R = 1 - 3 * ((beta * mu - log(beta + mu)) / (beta ** 2))**2 / 2
    S = 3 * log(beta + mu) / (2 * beta * mu) - 1.0_dp/2

    ex = ex * R
    Vx = Vx * S
end if
exc = ex + ec
Vxc = Vx + Vc

contains

    real(dp) function get_Y(y, b, c)
    real(dp), intent(in) :: y, b, c
    get_Y = y**2 + b*y + c
    end function

end subroutine

elemental subroutine xc_pz(n, exc, Vxc)
! Calculates XC LDA density and potential from the charge density "n".
! Uses the Perdew Zunger [1] parametrization.
!
! [1] Perdew, J. P., & Zunger, A. (1981). Self-interaction correction to
! density-functional approximations for many-electron systems. Physical Review
! B, 23(10), 5048–5079.
real(dp), intent(in) :: n ! charge density
real(dp), intent(out) :: exc ! XC density
real(dp), intent(out) :: Vxc ! XC potential

real(dp), parameter :: gam = -0.1423_dp
real(dp), parameter :: beta1 = 1.0529_dp
real(dp), parameter :: beta2 = 0.3334_dp
real(dp), parameter :: A =  0.0311_dp
real(dp), parameter :: B = -0.048_dp
real(dp), parameter :: C =  0.0020_dp
real(dp), parameter :: D = -0.0116_dp
real(dp) :: ex, ec, Vx, Vc, rs, sqrt_rs, log_rs

ex = -3/(4*pi) * (3*pi**2*n)**(1.0_dp/3)
Vx = 4*ex/3

rs = (3/(4*pi*n))**(1.0_dp/3)
if (rs >= 1) then
    sqrt_rs = sqrt(rs)
    ec = gam / (1+beta1*sqrt_rs+beta2*rs)
    Vc = ec * (1+7*beta1*sqrt_rs/6 + 4*beta2*rs/3) / &
        (1+beta1*sqrt_rs + beta2*rs)
else
    log_rs = log(rs)
    ec = A*log_rs + B + C*rs*log_rs + D*rs
    Vc = A*log_rs + (B-A/3) + 2*C*rs*log_rs/3 + (2*D-C)*rs/3
end if

exc = ex + ec
Vxc = Vx + Vc
end subroutine

end module
