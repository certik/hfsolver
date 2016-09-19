module ffte
use iso_fortran_env, only: real64
use types, only: dp
implicit none
private
public fft3_inplace, ifft3_inplace, pfft3, pfft3_init, pifft3

! FFTE is using REAL*8, which is 8 bytes real = 64 bit real = real64
integer, parameter:: dp_ffte=kind(1._real64)

interface

    subroutine FACTOR(N, IP)
    integer, intent(in) :: n
    integer, intent(out) :: ip(*)
    end subroutine

    subroutine ZFFT3D(A,NX,NY,NZ,IOPT)
    import :: dp_ffte
    complex(dp_ffte), intent(inout) :: A(*)
    integer, intent(in) :: NX, NY, NZ, IOPT
    end subroutine

    subroutine PZFFT3DV(A,B,NX,NY,NZ,ICOMMY,ICOMMZ,NPUY,NPUZ,IOPT)
    import :: dp_ffte
    complex(dp_ffte), intent(inout) :: A(*), B(*)
    integer, intent(in) :: NX, NY, NZ, ICOMMY, ICOMMZ, NPUY, NPUZ, IOPT
    end subroutine
end interface

contains

subroutine fft3_inplace(x)
complex(dp), intent(inout) :: x(:, :, :)
call zfft3d(x, size(x, 1), size(x, 2), size(x, 3), 0)
call zfft3d(x, size(x, 1), size(x, 2), size(x, 3), -1)
end subroutine

subroutine ifft3_inplace(x)
complex(dp), intent(inout) :: x(:, :, :)
call zfft3d(x, size(x, 1), size(x, 2), size(x, 3), 0)
call zfft3d(x, size(x, 1), size(x, 2), size(x, 3), 1)
x = x * product(shape(x))
end subroutine

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
integer, intent(in) :: commy, commz
integer, intent(in) :: Ng(:), nsub(:)
call pzfft3dv(x, y, Ng(1), Ng(2), Ng(3), commy, commz, nsub(2), nsub(3), -1)
end subroutine

subroutine pifft3(x, y, commy, commz, Ng, nsub)
complex(dp), intent(inout) :: x(:, :, :) ! input, will get destroyed
complex(dp), intent(out) :: y(:, :, :)   ! output, will get divided by N
integer, intent(in) :: commy, commz
integer, intent(in) :: Ng(:), nsub(:)
call pzfft3dv(x, y, Ng(1), Ng(2), Ng(3), commy, commz, nsub(2), nsub(3), 1)
end subroutine

end module
