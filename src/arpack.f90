module arpack
use lapack, only: daxpy, dscal, dnrm2
use types, only: dp
use utils, only: stop_error
implicit none
private
public eig, peig

interface
    subroutine mul_matrix_vector(x, y)
    import :: dp
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    end subroutine
end interface

! This is the precision that ARPACK "d" routines were compiled with.
integer, parameter:: adp=kind(0.d0)

interface

    subroutine dsaupd &
         ( ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, &
           ipntr, workd, workl, lworkl, info )
    import :: adp
    character bmat*1, which*2
    integer  ido, info, ldv, lworkl, n, ncv, nev
    real(adp) tol
    integer  iparam(11), ipntr(11)
    real(adp) resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
    end subroutine

    subroutine pdsaupd &
         ( comm, ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, iparam, &
           ipntr, workd, workl, lworkl, info )
    import :: adp
    integer :: comm
    character bmat*1, which*2
    integer  ido, info, ldv, lworkl, n, ncv, nev
    real(adp) tol
    integer  iparam(11), ipntr(11)
    real(adp) resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
    end subroutine

    subroutine dseupd (rvec  , howmny, select, d    , &
                       z     , ldz   , sigma , bmat , &
                       n     , which , nev   , tol  , &
                       resid , ncv   , v     , ldv  , &
                       iparam, ipntr , workd , workl, &
                       lworkl, info )
    import :: adp
    character  bmat*1, howmny*1, which*2
    logical    rvec
    integer    info, ldz, ldv, lworkl, n, ncv, nev
    real(adp) sigma, tol
    integer    iparam(7), ipntr(11)
    logical    select(ncv)
    real(adp) d(nev)     , resid(n)  , v(ldv,ncv), &
             z(ldz, nev), workd(2*n), workl(lworkl)
    end subroutine

    subroutine pdseupd (comm, rvec  , howmny, select, d    , &
                       z     , ldz   , sigma , bmat , &
                       n     , which , nev   , tol  , &
                       resid , ncv   , v     , ldv  , &
                       iparam, ipntr , workd , workl, &
                       lworkl, info )
    import :: adp
    integer :: comm
    character  bmat*1, howmny*1, which*2
    logical    rvec
    integer    info, ldz, ldv, lworkl, n, ncv, nev
    real(adp) sigma, tol
    integer    iparam(7), ipntr(11)
    logical    select(ncv)
    real(adp) d(nev)     , resid(n)  , v(ldv,ncv), &
             z(ldz, nev), workd(2*n), workl(lworkl)
    end subroutine


end interface

contains


subroutine eig(n, nev, ncv, which, matvec, d, v)
integer, intent(in) :: n ! Dimension of the matrix A
integer, intent(in) :: nev ! Number of eigenvalues to compute
integer, intent(in) :: ncv ! How many Lanczos vectors to generate
! Constraints: nev < ncv <= n
procedure(mul_matrix_vector) :: matvec ! returns y = A*x
character(2), intent(in) :: which
real(dp), intent(out) :: d(:), v(:,:)

real(dp) :: workl(ncv*(ncv+8)), workd(3*n), resid(n), tol, sigma
integer :: iparam(11), ipntr(11), ido, info, ierr, ishfts, maxitr, mode1, nconv
logical :: select(ncv), rvec
character(1) :: bmat
bmat  = 'I'
tol = 0
info = 0
ido = 0
ishfts = 1
maxitr = 300
mode1 = 1
iparam(1) = ishfts
iparam(3) = maxitr
iparam(7) = mode1
do
    call dsaupd ( ido, bmat, n, which, nev, tol, resid, &
                  ncv, v, size(v,1), iparam, ipntr, workd, workl, &
                  size(workl), info )
    if (ido == 99) exit
    if (ido == -1 .or. ido == 1) then
        call matvec(workd(ipntr(1):ipntr(1)+n-1), workd(ipntr(2):ipntr(2)+n-1))
    else
        call stop_error("Incorrect 'ido'.")
    end if
end do
if ( info .lt. 0 ) then
    print *, ' '
    print *, ' Error with _saupd, info = ', info
    print *, ' Check documentation in _saupd '
    print *, ' '
    call stop_error("")
end if
rvec = .true.
call dseupd ( rvec, 'All', select, d, v, size(v,1), sigma, &
    bmat, n, which, nev, tol, resid, ncv, v, size(v,1), &
    iparam, ipntr, workd, workl, size(workl), ierr )
if ( ierr .ne. 0) then
    print *, ' '
    print *, ' Error with _seupd, info = ', ierr
    print *, ' Check the documentation of _seupd. '
    print *, ' '
    call stop_error("")
end if
nconv =  iparam(5)
if (nconv < nev) &
    call stop_error("Some eigenvalues did not converge (nconv < nev).")
if ( info .eq. 1) then
    print *, ' '
    print *, ' Maximum number of iterations reached.'
    print *, ' '
    call stop_error("")
end if
if ( info .eq. 3) then
    print *, ' '
    print *, ' No shifts could be applied during implicit', &
             ' Arnoldi update, try increasing NCV.'
    print *, ' '
    call stop_error("")
end if
print *, ' '
print *, ' _SSIMP '
print *, ' ====== '
print *, ' '
print *, ' Size of the matrix is ', n
print *, ' The number of Ritz values requested is ', nev
print *, ' The number of Arnoldi vectors generated', ' (NCV) is ', ncv
print *, ' What portion of the spectrum: ', which
print *, ' The number of converged Ritz values is ', nconv
print *, ' The number of Implicit Arnoldi update iterations taken is ',iparam(3)
print *, ' The number of OP*x is ', iparam(9)
print *, ' The convergence criterion is ', tol
print *, ' '
end subroutine

subroutine peig(comm, myid, n, nev, ncv, which, matvec, d, v)
integer, intent(in) :: comm, myid ! MPI communicator and my ID
integer, intent(in) :: n ! Dimension of the matrix A
integer, intent(in) :: nev ! Number of eigenvalues to compute
integer, intent(in) :: ncv ! How many Lanczos vectors to generate
! Constraints: nev < ncv <= n
procedure(mul_matrix_vector) :: matvec ! returns y = A*x
character(2), intent(in) :: which
real(dp), intent(out) :: d(:), v(:,:)

real(dp) :: workl(ncv*(ncv+8)), workd(3*n), resid(n), tol, sigma
integer :: iparam(11), ipntr(11), ido, info, ierr, ishfts, maxitr, mode1, nconv
logical :: select(ncv), rvec
character(1) :: bmat
logical :: verbose
verbose = .false.
bmat  = 'I'
tol = 0
info = 0
ido = 0
ishfts = 1
maxitr = 300
mode1 = 1
iparam(1) = ishfts
iparam(3) = maxitr
iparam(7) = mode1
do
    call pdsaupd ( comm, ido, bmat, n, which, nev, tol, resid, &
                  ncv, v, size(v,1), iparam, ipntr, workd, workl, &
                  size(workl), info )
    if (ido == 99) exit
    if (ido == -1 .or. ido == 1) then
        call matvec(workd(ipntr(1):ipntr(1)+n-1), workd(ipntr(2):ipntr(2)+n-1))
    else
        call stop_error("Incorrect 'ido'.")
    end if
end do
if ( info .lt. 0 ) then
    print *, ' '
    print *, ' Error with _saupd, info = ', info
    print *, ' Check documentation in _saupd '
    print *, ' '
    call stop_error("")
end if
rvec = .true.
call pdseupd ( comm, rvec, 'All', select, d, v, size(v,1), sigma, &
    bmat, n, which, nev, tol, resid, ncv, v, size(v,1), &
    iparam, ipntr, workd, workl, size(workl), ierr )
if ( ierr .ne. 0) then
    print *, ' '
    print *, ' Error with _seupd, info = ', ierr
    print *, ' Check the documentation of _seupd. '
    print *, ' '
    call stop_error("")
end if
nconv =  iparam(5)
if (nconv < nev) &
    call stop_error("Some eigenvalues did not converge (nconv < nev).")
if ( info .eq. 1) then
    print *, ' '
    print *, ' Maximum number of iterations reached.'
    print *, ' '
    call stop_error("")
end if
if ( info .eq. 3) then
    print *, ' '
    print *, ' No shifts could be applied during implicit', &
             ' Arnoldi update, try increasing NCV.'
    print *, ' '
    call stop_error("")
end if
if (myid == 0 .and. verbose) then
    print *, ' '
    print *, ' _SSIMP '
    print *, ' ====== '
    print *, ' '
    print *, ' The number of Ritz values requested is ', nev
    print *, ' The number of Arnoldi vectors generated', ' (NCV) is ', ncv
    print *, ' What portion of the spectrum: ', which
    print *, ' The number of converged Ritz values is ', nconv
    print *, ' The number of Implicit Arnoldi update iterations taken is ',iparam(3)
    print *, ' The number of OP*x is ', iparam(9)
    print *, ' The convergence criterion is ', tol
    print *, ' '
end if
end subroutine

end module
