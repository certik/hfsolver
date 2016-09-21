module pffte
use types, only: dp
use ffte, only: dp_ffte, factor
use utils, only: stop_error
use mpi2, only: mpi_barrier
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

subroutine pfft3_init(myid, comm_all, Ng, nsub)
integer, intent(in) :: myid, comm_all, Ng(:), nsub(:)
complex(dp) :: x(1), y(1)
integer :: commy, commz, ierr
integer :: Ng_local(3), F(3)
if (myid==0) then
    if (.not. all(nsub > 0)) then
        call stop_error("nsub must be positive")
    end if
    Ng_local = Ng / nsub
    if (nsub(1) /= 1) then
        call stop_error("nsub(1) must be equal to 1")
    end if
    if (.not. all(Ng_local * nsub == Ng)) then
        call stop_error("Ng must be equal to Ng_local * nsub")
    end if
    if (mod(Ng(1), nsub(2)) + mod(Ng(1), nsub(3)) + mod(Ng(2), nsub(3)) &
        + mod(Ng(3), nsub(2)) > 0) then
        call stop_error("Ng must be divisible by any permutation of nsub")
    end if
    call factor(Ng(1), F)
    if (Ng(1) /= 2**F(1)*3**F(2)*5**F(3)) then
        call stop_error("Ng(1) must be 2^a * 3^b * 5^c;  a,b,c=0,1,2,3,...")
    end if
    call factor(Ng(2), F)
    if (Ng(2) /= 2**F(1)*3**F(2)*5**F(3)) then
        call stop_error("Ng(2) must be 2^a * 3^b * 5^c;  a,b,c=0,1,2,3,...")
    end if
    call factor(Ng(3), F)
    if (Ng(3) /= 2**F(1)*3**F(2)*5**F(3)) then
        call stop_error("Ng(3) must be 2^a * 3^b * 5^c;  a,b,c=0,1,2,3,...")
    end if
end if
call mpi_barrier(comm_all, ierr)
! Only Ng is relevant if iopt == 0
! pzfft3dv divides by nsub inside, but does not store the result (the above
! checks make sure we do not divide by 0)
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
