module pffte
use types, only: dp
use ffte, only: dp_ffte
implicit none
private
public pfft3_init, pfft3, pifft3

interface
    subroutine PZFFT3DV(A,B,NX,NY,NZ,ICOMMY,ICOMMZ,NPUY,NPUZ,IOPT)
    import :: dp_ffte
    complex(dp_ffte), intent(inout) :: A(*), B(*)
    integer, intent(in) :: NX, NY, NZ, ICOMMY, ICOMMZ, NPUY, NPUZ, IOPT
    end subroutine
end interface

contains

subroutine pfft3_init(Ng)
integer, intent(in) :: Ng(:)
complex(dp) :: x(1), y(1)
integer :: commy, commz, nsub(3)
! Only Ng is relevant if iopt == 0
nsub = 1 ! pzfft3dv divides by nsub inside, but does not store the result
call pzfft3dv(x, y, Ng(1), Ng(2), Ng(3), commy, commz, nsub(2), nsub(3), 0)
end subroutine

subroutine pfft3(x, y, commy, commz, Ng, nsub)
complex(dp), intent(inout) :: x(:, :, :) ! input, will get destroyed
complex(dp), intent(out) :: y(:, :, :)   ! output
integer, intent(in) :: commy, commz ! communicators in y, z directions
integer, intent(in) :: Ng(:) ! Total (global) number of PW in each direction
integer, intent(in) :: nsub(:) ! Number of subdomains in each direction
call pzfft3dv(x, y, Ng(1), Ng(2), Ng(3), commy, commz, nsub(2), nsub(3), -1)
end subroutine

subroutine pifft3(x, y, commy, commz, Ng, nsub)
complex(dp), intent(inout) :: x(:, :, :) ! input, will get destroyed
complex(dp), intent(out) :: y(:, :, :)   ! output, will get divided by N
integer, intent(in) :: commy, commz ! communicators in y, z directions
integer, intent(in) :: Ng(:) ! Total (global) number of PW in each direction
integer, intent(in) :: nsub(:) ! Number of subdomains in each direction
call pzfft3dv(x, y, Ng(1), Ng(2), Ng(3), commy, commz, nsub(2), nsub(3), 1)
! pzfft3dv divides by N, so we multply by N here to cancel that
y = y * product(Ng)
end subroutine

end module