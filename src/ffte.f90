module ffte
use iso_fortran_env, only: real64
use types, only: dp
implicit none
private
public fft3_inplace, ifft3_inplace

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

end module
