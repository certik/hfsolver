module xc

! Exchange and correlation potentials

use types, only: dp
use constants, only: pi
implicit none
private
public get_Vxc

contains

subroutine get_Vxc(R, rho, relat, c, exc, Vxc)
real(dp), intent(in) :: R(:) ! radial grid
real(dp), intent(in) :: rho(:) ! charge density
logical, intent(in) :: relat ! .true. return RLDA, otherwise LDA
real(dp), intent(in) :: c ! speed of light
real(dp), intent(out) :: Vxc(:), exc(:)

integer :: i
do i = 1, size(R)
    call getvxc_scalar(rho(i), relat, c, exc(i), Vxc(i))
end do
end subroutine

real(dp) function get_Y(y, b, c)
real(dp), intent(in) :: y, b, c
get_Y = y**2 + b*y + c
end function

subroutine getvxc_scalar(n, relat, c_light, exc, Vxc)
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

if (n == 0) then
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
end subroutine

end module
