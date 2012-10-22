module linalg_feast
use types, only: dp
use utils, only: stop_error, assert
use feast, only: feastinit, dfeast_syev, dfeast_srci
use lapack, only: xerbla, zcopy, zgetrf, zgetrs, dgemm, dsymm
use petsc_, only: petsc_init, petsc_finalize, solve
implicit none
private
public eigh

interface eigh
    module procedure deigh_generalized
    module procedure deigh_simple
end interface eigh


contains

subroutine deigh_generalized(A, B, Emin, Emax, M0, lam, c)
real(dp), intent(in) :: A(:, :), B(:, :), Emin, Emax
integer, intent(in) :: M0
real(dp), allocatable, intent(out) :: lam(:), c(:, :)
integer :: feastparam(64)
real(dp), allocatable :: res(:), lam_(:), c_(:, :)
real(dp) :: epsout
integer :: i, loop, info, M, N, LDA
N = size(A, 1)
LDA = N
allocate(lam_(M0), c_(N, M0), res(M0))

call feastinit(feastparam)
feastparam(1)=1 !! change from default value
feastparam(2) = 4  !! Nq
feastparam(3) = 5  !! accuracy: 1e-x
call dfeast_sygv('L',N,A,LDA,B,LDA,feastparam,epsout,loop,Emin,Emax,M0, &
    lam_,c_,M,res,info)

if (info /= 0) then
    print *,'FEAST OUTPUT INFO', info
    call stop_error("info /= 0")
end if

print *, '# Search interval [Emin,Emax]', Emin, Emax
print *, '# mode found/subspace', M, M0
print *, '# iterations', loop
print *, 'TRACE', sum(lam_(1:M))
print *, 'Relative error on the Trace', epsout
print *, 'Eigenvalues/Residuals'
do i = 1, M
    print *, i, lam_(i), res(i)
end do

allocate(lam(M), c(N, M))
lam = lam_(:M)
c = c_(:, :M)

end subroutine


subroutine deigh_simple(A, Emin, Emax, M0, lam, c)
real(dp), intent(in) :: A(:, :), Emin, Emax
integer, intent(in) :: M0
real(dp), intent(out) :: lam(:), c(:, :)
integer :: feastparam(64)
real(dp), allocatable :: res(:)
real(dp) :: epsout
integer :: i, loop, info, M, N, LDA
N = size(A, 1)
LDA = N
allocate(res(M0))

call feastinit(feastparam)
feastparam(1)=1 !! change from default value
call dfeast_syev('F',N,A,LDA,feastparam,epsout,loop,Emin,Emax,M0,lam,c, &
    M,res,info)

if (info /= 0) then
    print *,'FEAST OUTPUT INFO', info
    call stop_error("info /= 0")
end if

print *, '# Search interval [Emin,Emax]', Emin, Emax
print *, '# mode found/subspace', M, M0
print *, '# iterations', loop
print *, 'TRACE', sum(lam(1:M))
print *, 'Relative error on the Trace', epsout
print *, 'Eigenvalues/Residuals'
do i = 1, M
    print *, i, lam(i), res(i)
end do

if (M < M0) call stop_error("M < M0")

end subroutine

subroutine dfeast_sygv(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        REAL DOUBLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) REAL DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                        if option fpm(6)=1 selected                           
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2009-2012
  ! ====================================================================
  character(len=1) :: UPLO
  integer :: N,LDA,LDB
  double precision,dimension(LDA,*):: A
  double precision,dimension(LDB,*):: B
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  integer :: M0
  double precision,dimension(*)  :: E
  double precision,dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
  integer :: ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::workc,Az
  double precision, dimension(:,:),pointer ::work,Aq,Sq
  integer, dimension(:),pointer ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)

  INFO = 0
  IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
     INFO=-101
  ELSE IF ( N<=0 ) THEN
     INFO = -102
  ELSE IF(LDA<N ) THEN
     INFO = -104
  ELSE IF(LDB<N ) THEN
     INFO = -106
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DFEAST_SYGV', -INFO+100 )
     RETURN
  END IF
  infoloc=0
  allocate(Aq(M0,M0), Sq(M0, M0), work(N, M0), workc(N, M0), ipivloc(N))
  allocate(Az(N, N))
  ierr = petsc_init()

  if (infoloc/=0) then
     info=-1
     return
  end if
!!$
  ijob=-1 ! initialization 
  do while (ijob/=0)
      call dfeast_srci(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
     select case(ijob)
     case(10) !! factorize (zeB-A)
        Az(1:N,1:N)=Ze*B(1:N,1:N)-A(1:N,1:N)*ONEc

        print *, "Factorization (zgetrf) N =", N
        call ZGETRF(N,N,Az,N,IPIVloc,INFOloc)
        print *, "Done"
        if (infoloc/=0) then
           info=-2
           return
        end if



     case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:M0) result in to workc


        print *, "Solve (zgetrs) N =", N
        call ZGETRS( 'N', N, M0, Az, N, IPIVloc, workc, N, INFOloc )
        print *, "Done"
        if (infoloc/=0) then
           info=-2
           return
        end if


     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if ((UPLO=='F').or.(UPLO=='f')) then
           call DGEMM('N','N',N,fpm(25),N,DONE,A,LDA,X(1,fpm(24)),N,DZERO,work(1:N,fpm(24)),N)
        else
           call DSYMM ('L', UPLO, N, fpm(25), DONE, A, LDA, X(1,fpm(24)), N, DZERO,work(1:N,fpm(24)), N)
        endif

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if ((UPLO=='F').or.(UPLO=='f')) then
           call DGEMM('N','N',N,fpm(25),N,DONE,B,LDB,X(1,fpm(24)),N,DZERO,work(1:N,fpm(24)),N)
        else
           call DSYMM ('L', UPLO, N, fpm(25), DONE, B, LDB, X(1,fpm(24)), N, DZERO,work(1:N,fpm(24)), N)
        endif

     end select
  end do

  ierr = petsc_finalize()


end subroutine dfeast_sygv

end module
