module pksdft_fft
use constants, only: pi
use types, only: dp
use utils, only: assert, stop_error, clock, allocate_mold
use pofdft_fft, only: preal2fourier, pfourier2real
use arpack, only: peig
implicit none
private
public solve_schroedinger

contains

subroutine solve_schroedinger(myid, comm_all, commy, commz, Ng, nsub, Vloc, &
        L, G2, cutfn, nev, ncv, eigs, orbitals)
integer, intent(in) :: myid, comm_all, commy, commz, Ng(3), nsub(3), nev, ncv
real(dp), intent(in) :: Vloc(:,:,:) ! Local effective potential
real(dp), intent(in) :: G2(:,:,:), L(:), cutfn(:,:,:)
real(dp), intent(out) :: eigs(:) ! eigs(nev)
! orbitals(Ng_local(1),Ng_local(2),Ng_local(3),nev)
real(dp), intent(out) :: orbitals(:,:,:,:)

real(dp), allocatable :: d(:), v(:,:)
integer :: Ng_local(3), n, i
real(dp) :: t1, t2
logical :: verbose
verbose = .false.
Ng_local = Ng / nsub
n = product(Ng_local)
allocate(v(n,ncv), d(ncv))
call cpu_time(t1)
call peig(comm_all, myid, n, nev, ncv, "SA", av, d, v)
call cpu_time(t2)
eigs = d(:nev)+10
orbitals = reshape(v(:,:nev), [Ng_local(1),Ng_local(2),Ng_local(3),nev]) &
    * sqrt(product(Ng/L))
if (myid == 0 .and. verbose) then
    print *, "Arpack Time:", t2-t1
end if

contains

    subroutine av(x, y)
    ! Compute y = A*x
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    complex(dp), dimension(Ng_local(1),Ng_local(2),Ng_local(3)) :: &
        psi, psiG, psiG_vloc
    call preal2fourier(reshape(x, [Ng_local(1),Ng_local(2),Ng_local(3)]), &
        psiG, commy, commz, Ng, nsub)
    psiG = psiG * cutfn
    ! psiG is our starting point

    ! Apply kinetic and potential
    call pfourier2real(psiG, psi, commy, commz, Ng, nsub)
    call preal2fourier((Vloc-10)*psi, psiG_vloc, commy, commz, Ng, nsub)
    psiG_vloc = psiG_vloc * cutfn

    psiG = G2/2*psiG + psiG_vloc

    ! Convert to real space at the end
    call pfourier2real(psiG, psi, commy, commz, Ng, nsub)
    y = reshape(real(psi,dp), [product(Ng_local)])
    end

end subroutine

end module
