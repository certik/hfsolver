module ofdft

! Routines for orbital free density functional theory. These routines are
! independent of any particular discretization.

use types, only: dp
use constants, only: pi
implicit none
private
public f, read_pseudo

contains

real(dp) elemental function f(y, deriv)
! Function f(y) from Appendix A in [1].
!
! [1] Perrot, F. (1979). Gradient correction to the statistical electronic free
! energy at nonzero temperatures: Application to equation-of-state
! calculations. Physical Review A, 20(2), 586–594.
real(dp), intent(in) :: y ! must be positive
! if deriv == .true. compute df/dy instead. Default .false.
logical, intent(in), optional :: deriv
real(dp), parameter :: y0 = 3*pi/(4*sqrt(2._dp))
real(dp), parameter :: c(8) = [-0.8791880215_dp, 0.1989718742_dp, &
    0.1068697043e-2_dp, -0.8812685726e-2_dp, 0.1272183027e-1_dp, &
    -0.9772758583e-2_dp, 0.3820630477e-2_dp, -0.5971217041e-3_dp]
real(dp), parameter :: d(9) = [0.7862224183_dp, -0.1882979454e1_dp, &
    0.5321952681_dp, 0.2304457955e1_dp, -0.1614280772e2_dp, &
    0.5228431386e2_dp, -0.9592645619e2_dp, 0.9462230172e2_dp, &
    -0.3893753937e2_dp]
real(dp) :: u
integer :: i
logical :: deriv_
deriv_ = .false.
if (present(deriv)) deriv_ = deriv

if (.not. deriv_) then
    if (y <= y0) then
        f = log(y)
        do i = 0, 7
            f = f + c(i+1) * y**i
        end do
    else
        u = y**(2._dp / 3)
        f = 0
        do i = 0, 8
            f = f + d(i+1) / u**(2*i-1)
        end do
        ! Note: Few terms in [1] have "y" instead of "u" in them for y > y0, but
        ! that is obviously a typo.
    end if
else
    if (y <= y0) then
        f = 1 / y
        do i = 0, 6
            f = f + (i+1) * c(i+2) * y**i
        end do
    else
        u = y**(2._dp / 3)
        f = 0
        do i = 0, 8
            f = f - (2*i-1) * d(i+1) / u**(2*i)
        end do
        f = f * 2._dp/3 / y**(1._dp/3)
    end if
end if
end function

subroutine read_pseudo(filename, R, V, Z, Ediff)
! Reads the pseudopotential from the file 'filename'.
character(len=*), intent(in) :: filename   ! File to read from, e.g. "H.pseudo"
real(dp), allocatable, intent(out) :: R(:) ! radial grid [0, Rcut]
! potential on the radial grid. The potential smoothly changes into -1/R for
! r > Rcut, where Rcut = R(size(R)) is the cut-off radius
real(dp), allocatable, intent(out) :: V(:)
real(dp), intent(out) :: Z     ! Nuclear charge
real(dp), intent(out) :: Ediff ! The energy correction
real(dp) :: Rcut
integer :: N, i, u
open(newunit=u, file=filename, status="old")
read(u, *) Z, N, Rcut, Ediff
allocate(R(N), V(N))
do i = 1, N
    read(u, *) R(i), V(i)
end do
close(u)
R = R * Rcut
end subroutine

end module
