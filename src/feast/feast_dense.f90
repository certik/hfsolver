!=========================================================================================
!Copyright (c) 2009-2012, The Regents of the University of Massachusetts, Amherst.
!Developed by E. Polizzi
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without modification, 
!are permitted provided that the following conditions are met:
!
!1. Redistributions of source code must retain the above copyright notice, this list of conditions 
!   and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
!   and the following disclaimer in the documentation and/or other materials provided with the distribution.
!3. Neither the name of the University nor the names of its contributors may be used to endorse or promote
!    products derived from this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
!BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
!ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
!EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
!LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
!IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!==========================================================================================



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! FEAST PREDEFINED DENSE INTERFACES !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! List of routines:
!-------------------

!{S,D,C,Z}FEAST_{SY,HE}{EV,GV}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Symmetric eigenvalue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! single/double precision
! real symmetric 
!{s,d}feast_sygv ! generalized
!{s,d}feast_syev ! standard

! complex Hermitian
!{c,z}feast_hegv ! generalized
!{c,z}feast_heev ! standard
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$

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
  implicit none
  include 'f90_noruntime_interface.fi'
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc,i,s
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
  call wallocate_2d(Aq,M0,M0,infoloc)
  call wallocate_2d(Sq,M0,M0,infoloc)
  call wallocate_2d(work,N,M0,infoloc)
  call wallocate_2z(workc,N,M0,infoloc)
  call wallocate_1i(ipivloc,N,infoloc)
  call wallocate_2z(Az,N,N,infoloc)

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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION

        if ((UPLO=='L').or.(UPLO=='l')) then
           do i=1,N-1
              s=N-i
              Az(i:N,i)=Ze*B(i:N,i)-A(i:N,i)*ONEC
              call ZCOPY(s,Az(i+1,i),1,Az(i,i+1),N)
           enddo
           Az(N,N)=Ze*B(N,N)-A(N,N)*ONEC

        elseif ((UPLO=='U').or.(UPLO=='u')) then
           do i=1,N-1
              s=N-i
              Az(i,i:N)=Ze*B(i,i:N)-A(i,i:N)*ONEC
              call ZCOPY(s,Az(i,i+1),N,Az(i+1,i),1)
           enddo
           Az(N,N)=Ze*B(N,N)-A(N,N)*ONEC

        else ! full 
           Az(1:N,1:N)=Ze*B(1:N,1:N)-A(1:N,1:N)*ONEc

        end if



        call ZGETRF(N,N,Az,N,IPIVloc,INFOloc)     
        if (infoloc/=0) then
           info=-2
           return
        end if



     case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:M0) result in to workc


        call ZGETRS( 'N', N, M0, Az, N, IPIVloc, workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if


     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if ((UPLO=='F').or.(UPLO=='f')) then
           call DGEMM('N','N',N,fpm(25),N,DONE,A,LDA,X(1,fpm(24)),N,DZERO,work(1,fpm(24)),N)
        else
           call DSYMM ('L', UPLO, N, fpm(25), DONE, A, LDA, X(1,fpm(24)), N, DZERO,work(1,fpm(24)), N)
        endif

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if ((UPLO=='F').or.(UPLO=='f')) then
           call DGEMM('N','N',N,fpm(25),N,DONE,B,LDB,X(1,fpm(24)),N,DZERO,work(1,fpm(24)),N)
        else
           call DSYMM ('L', UPLO, N, fpm(25), DONE, B, LDB, X(1,fpm(24)), N, DZERO,work(1,fpm(24)), N)
        endif

     end select
  end do



  call wdeallocate_2d(Aq)
  call wdeallocate_2d(Sq)
  call wdeallocate_2d(work)
  call wdeallocate_2z(workc)
  call wdeallocate_2z(Az)
  call wdeallocate_1i(ipivloc)


end subroutine dfeast_sygv



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine dfeast_syev(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A REAL SYMMETRIC:: DENSE FORMAT 
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
  implicit none
  include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,LDA
  double precision,dimension(LDA,*):: A
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc,i,s
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::workc,Az
  double precision, dimension(:,:),pointer ::work,Aq,Sq
  integer, dimension(:),pointer ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)

  INFO = 0
  IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
     INFO=-101
  ELSE IF( N<=0 ) THEN
     INFO = -102
  ELSE IF(LDA<N ) THEN
     INFO = -104
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'DFEAST_SYEV', -INFO+100 )
     RETURN
  END IF

  infoloc=0
  call wallocate_2d(Aq,M0,M0,infoloc)
  call wallocate_2d(Sq,M0,M0,infoloc)
  call wallocate_2d(work,N,M0,infoloc)
  call wallocate_2z(workc,N,M0,infoloc)
  call wallocate_2z(Az,N,N,infoloc)
  call wallocate_1i(ipivloc,N,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if


  ijob=-1 ! initialization

  do while (ijob/=0)
     call dfeast_srci(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)

     select case(ijob)
     case(10) !! factorize (zeB-A)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION

        if ((UPLO=='L').or.(UPLO=='l')) then
           do i=1,N-1
              s=N-i
              Az(i:N,i)=-A(i:N,i)*ONEC
              Az(i,i)=Az(i,i)+Ze
              call ZCOPY(s,Az(i+1,i),1,Az(i,i+1),N)
           enddo
           Az(N,N)=Ze-A(N,N)*ONEC

        elseif ((UPLO=='U').or.(UPLO=='u')) then
           do i=1,N-1
              s=N-i
              Az(i,i:N)=-A(i,i:N)*ONEC
              Az(i,i)=Az(i,i)+Ze
              call ZCOPY(s,Az(i,i+1),N,Az(i+1,i),1)
           enddo
           Az(N,N)=Ze-A(N,N)*ONEC


        else

           Az(1:N,1:N)=-A(1:N,1:N)*ONEC
           do i=1,N
              Az(i,i)=Az(i,i)+Ze
           enddo

        end if


        call ZGETRF(N,N,Az,N,IPIVloc,INFOloc)     
        if (infoloc/=0) then
           info=-2
           return
        end if


     case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:M0) result in to workc

        call ZGETRS( 'N', N, M0, Az, N, IPIVloc, workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if


     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if ((UPLO=='F').or.(UPLO=='f')) then
           call DGEMM('N','N',N,fpm(25),N,DONE,A,LDA,X(1,fpm(24)),N,DZERO,work(1,fpm(24)),N)
        else
           call DSYMM ('L', UPLO, N, fpm(25), DONE, A, LDA, X(1,fpm(24)), N, DZERO,work(1,fpm(24)), N)
        endif

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        call DLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )

     end select
  end do


  call wdeallocate_2d(Aq)
  call wdeallocate_2d(Sq)
  call wdeallocate_2d(work)
  call wdeallocate_2z(workc)
  call wdeallocate_2z(Az)
  call wdeallocate_1i(ipivloc)


end subroutine dfeast_syev


!!$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine zfeast_hegv(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B 
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
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                        if option fpm(6)=1 selected                           
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2009-2012
  ! ====================================================================
  implicit none
  include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,LDA,LDB
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
  complex(kind=(kind(1.0d0))),dimension(LDB,*):: B
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  integer :: M0
  double precision,dimension(*)  :: E
  complex(kind=(kind(1.0d0))),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc,i,s
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::work1,work2,Az,zAq,zSq
  complex(kind=(kind(1.0d0))), dimension(:),pointer ::ztmp
  integer, dimension(:),pointer ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)



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
     CALL XERBLA( 'ZFEAST_HEGV', -INFO+100 )
     RETURN
  END IF


  infoloc=0
  call wallocate_2z(zAq,M0,M0,infoloc)
  call wallocate_2z(zSq,M0,M0,infoloc)
  call wallocate_2z(work1,N,M0,infoloc)
  call wallocate_2z(work2,N,M0,infoloc)
  call wallocate_2z(Az,N,N,infoloc)
  call wallocate_1i(ipivloc,N,infoloc)
  call wallocate_1z(ztmp,N,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if


  fpm(11)=1 ! half-contour only- requires solving (zB-A)^H x=f
  !fpm(11)=2 ! full contour  (2 half-contour)  - does not require solving for the transpose-conjg


  ijob=-1 ! initialization


  do while (ijob/=0) 
     call zfeast_hrci(ijob,N,Ze,work1,work2,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
     select case(ijob)

     case(10) !! factorize (zeB-A)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION

        if ((UPLO=='L').or.(UPLO=='l')) then
           do i=1,N-1
              s=N-i
              Az(i:N,i)=Ze*B(i:N,i)
              call ZCOPY(s,Az(i+1,i),1,Az(i,i+1),N)
              call ZLACGV(s, Az(i,i+1), N )
              call ZSCAL(s,Ze/conjg(Ze),Az(i,i+1),N)
              Az(i:N,i)=Az(i:N,i)-A(i:N,i)
           enddo
           Az(N,N)=Ze*B(N,N)-A(N,N)*ONEC
           do i=1,N-1
              s=N-i
              call ZCOPY(s,A(i+1,i),1,ztmp(1),1)
              call ZLACGV(s,ztmp(1),1)
              call ZAXPY(s,-ONEC,ztmp(1),1,Az(i,i+1),N)
           enddo


        elseif ((UPLO=='U').or.(UPLO=='u')) then

           do i=1,N-1
              s=N-i
              Az(i,i:N)=Ze*B(i,i:N)
              call ZCOPY(s,Az(i,i+1),N,Az(i+1,i),1)
              call ZLACGV(s, Az(i+1,i), 1 )
              call ZSCAL(s,Ze/conjg(Ze),Az(i+1,i),1)
              Az(i,i:N)=Az(i,i:N)-A(i,i:N)
           enddo
           Az(N,N)=Ze*B(N,N)-A(N,N)*ONEC
           do i=1,N-1
              s=N-i
              call ZCOPY(s,A(i,i+1),N,ztmp(1),1)
              call ZLACGV(s,ztmp(1),1)
              call ZAXPY(s,-ONEC,ztmp(1),1,Az(i+1,i),1)
           enddo

        else ! full 
           Az(1:N,1:N)=Ze*B(1:N,1:N)-A(1:N,1:N)
        end if


        call ZGETRF(N,N,Az,N,IPIVloc,INFOloc)     
        if (infoloc/=0) then
           info=-2
           return
        end if

     case(11) !!solve the linear system (ZeB-A)x=work2(1:N,1:M0) result in to work2


        call ZGETRS( 'N', N, M0, Az, N, IPIVloc, work2, N, INFOloc )


     case(21) !!solve the linear system (ZeB-A)^H x=work2(1:N,1:M0) result in to work2

        call ZGETRS( 'C', N, M0, Az, N, IPIVloc, work2, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if

     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if ((UPLO=='F').or.(UPLO=='f')) then
           call ZGEMM('N','N',N,fpm(25),N,ONEC,A,LDA,X(1,fpm(24)),N,ZEROC,work1(1,fpm(24)),N)
        else
           call ZHEMM ('L', UPLO, N, fpm(25), ONEC, A, LDA, X(1,fpm(24)), N, ZEROC,work1(1,fpm(24)), N)
        endif

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if ((UPLO=='F').or.(UPLO=='f')) then
           call ZGEMM('N','N',N,fpm(25),N,ONEC,B,LDB,X(1,fpm(24)),N,ZEROC,work1(1,fpm(24)),N)
        else
           call ZHEMM ('L', UPLO, N, fpm(25), ONEC, B, LDB, X(1,fpm(24)), N, ZEROC,work1(1,fpm(24)), N)
        endif

     end select
  end do

  call wdeallocate_2z(zAq)
  call wdeallocate_2z(zSq)
  call wdeallocate_2z(work1)
  call wdeallocate_2z(work2)
  call wdeallocate_2z(Az)
  call wdeallocate_1i(ipivloc)
  call wdeallocate_1z(ztmp)


end subroutine zfeast_hegv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine zfeast_heev(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN :: DENSE FORMAT 
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL DOUBLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL DOUBLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       REAL DOUBLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) COMPLEX DOUBLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL DOUBLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                        if option fpm(6)=1 selected                           
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2009-2012
  ! ====================================================================
  implicit none
  include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,LDA
  complex(kind=(kind(1.0d0))),dimension(LDA,*):: A
  integer,dimension(*) :: fpm
  double precision :: epsout 
  integer :: loop
  double precision :: Emin,Emax
  integer :: M0
  double precision,dimension(*)  :: E
  complex(kind=(kind(1.0d0))),dimension(N,*):: X
  integer :: mode
  double precision,dimension(*)    :: res
  integer :: info
!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc,i,s
  complex(kind=(kind(1.0d0))) :: Ze
  complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::work1,work2,Az,zAq,zSq
  integer, dimension(:),pointer ::ipivloc
  double precision,parameter :: DONE=1.0d0,DZERO=0.0d0
  complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)


  INFO = 0
  IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
     INFO=-101
  ELSE IF( N<=0 ) THEN
     INFO = -102
  ELSE IF(LDA<N ) THEN
     INFO = -104
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'ZFEAST_HEEV', -INFO+100 )
     RETURN
  END IF

  infoloc=0
  call wallocate_2z(zAq,M0,M0,infoloc)
  call wallocate_2z(zSq,M0,M0,infoloc)
  call wallocate_2z(work1,N,M0,infoloc)
  call wallocate_2z(work2,N,M0,infoloc)
  call wallocate_2z(Az,N,N,infoloc)
  call wallocate_1i(ipivloc,N,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if

  fpm(11)=1 ! half-contour only- requires solving (zB-A)^H x=f
  !fpm(11)=2 ! full contour  (2 half-contour)  - does not require solving for the transpose-conjg


  ijob=-1 ! initialization


  do while (ijob/=0) 
     call zfeast_hrci(ijob,N,Ze,work1,work2,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
     select case(ijob)
     case(10) !! factorize (zeB-A)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION

        if ((UPLO=='L').or.(UPLO=='l')) then
           do i=1,N-1
              s=N-i
              Az(i:N,i)=-A(i:N,i)
              Az(i,i)=Az(i,i)+Ze
              call ZCOPY(s,Az(i+1,i),1,Az(i,i+1),N)
              call ZLACGV(s, Az(i,i+1), N )
           enddo
           Az(N,N)=Ze-A(N,N)*ONEC

        elseif ((UPLO=='U').or.(UPLO=='u')) then
           do i=1,N-1
              s=N-i
              Az(i,i:N)=-A(i,i:N)
              Az(i,i)=Az(i,i)+Ze
              call ZCOPY(s,Az(i,i+1),N,Az(i+1,i),1)
              call ZLACGV(s, Az(i+1,i), 1 )
           enddo
           Az(N,N)=Ze-A(N,N)*ONEC

        else

           Az(1:N,1:N)=-A(1:N,1:N)
           do i=1,N
              Az(i,i)=Az(i,i)+Ze
           enddo

        end if



        call ZGETRF(N,N,Az,N,IPIVloc,INFOloc)     
        if (infoloc/=0) then
           info=-2
           return
        end if


     case(11) !!solve the linear system (ZeB-A)x=work2(1:N,1:M0) result in to work2


        call ZGETRS( 'N', N, M0, Az, N, IPIVloc, work2, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if

     case(21) !!solve the linear system (ZeB-A)^H x=work2(1:N,1:M0) result in to work2

        call ZGETRS( 'C', N, M0, Az, N, IPIVloc, work2, N, INFOloc )

     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if ((UPLO=='F').or.(UPLO=='f')) then
           call ZGEMM('N','N',N,fpm(25),N,ONEC,A,LDA,X(1,fpm(24)),N,ZEROC,work1(1,fpm(24)),N)
        else
           call ZHEMM ('L', UPLO, N, fpm(25), ONEC, A, LDA, X(1,fpm(24)), N, ZEROC,work1(1,fpm(24)), N)
        endif

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work1(1,fpm(24)), N )

     end select
  end do

  call wdeallocate_2z(zAq)
  call wdeallocate_2z(zSq)
  call wdeallocate_2z(work1)
  call wdeallocate_2z(work2)
  call wdeallocate_2z(Az)
  call wdeallocate_1i(ipivloc)


end subroutine zfeast_heev


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine sfeast_sygv(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        REAL SINGLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) REAL SINGLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                        if option fpm(6)=1 selected                           
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2009-2012
  ! ====================================================================
  implicit none
  include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,LDA,LDB
  real,dimension(LDA,*):: A
  real,dimension(LDB,*):: B
  integer,dimension(*) :: fpm
  real :: epsout 
  integer :: loop
  real :: Emin,Emax
  integer :: M0
  real,dimension(*)  :: E
  real,dimension(N,*):: X
  integer :: mode
  real,dimension(*)    :: res
  integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc,i,s
  complex  :: Ze
  complex , dimension(:,:),pointer ::workc,Az
  real, dimension(:,:),pointer ::work,Aq,Sq
  integer, dimension(:),pointer ::ipivloc
  real,parameter :: SONE=1.0E0,SZERO=0.0E0
  complex ,parameter :: ONEC=(SONE,SZERO)

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
     CALL XERBLA( 'SFEAST_SYGV', -INFO+100 )
     RETURN
  END IF

  infoloc=0
  call wallocate_2s(Aq,M0,M0,infoloc)
  call wallocate_2s(Sq,M0,M0,infoloc)
  call wallocate_2s(work,N,M0,infoloc)
  call wallocate_2c(workc,N,M0,infoloc)
  call wallocate_2c(Az,N,N,infoloc)
  call wallocate_1i(ipivloc,N,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if

  ijob=-1 ! initialization
  do while (ijob/=0)
     call sfeast_srci(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)



     select case(ijob)
     case(10) !! factorize (zeB-A)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION

        if ((UPLO=='L').or.(UPLO=='l')) then
           do i=1,N-1
              s=N-i
              Az(i:N,i)=Ze*B(i:N,i)-A(i:N,i)*ONEC
              call CCOPY(s,Az(i+1,i),1,Az(i,i+1),N)
           enddo
           Az(N,N)=Ze*B(N,N)-A(N,N)*ONEC

        elseif ((UPLO=='U').or.(UPLO=='u')) then
           do i=1,N-1
              s=N-i
              Az(i,i:N)=Ze*B(i,i:N)-A(i,i:N)*ONEC
              call CCOPY(s,Az(i,i+1),N,Az(i+1,i),1)
           enddo
           Az(N,N)=Ze*B(N,N)-A(N,N)*ONEC

        else ! full 
           Az(1:N,1:N)=Ze*B(1:N,1:N)-A(1:N,1:N)*ONEc
        end if



        call CGETRF(N,N,Az,N,IPIVloc,INFOloc)     
        if (infoloc/=0) then
           info=-2
           return
        end if


     case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:M0) result in to workc

        call CGETRS( 'N', N, M0, Az, N, IPIVloc, workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if


     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if ((UPLO=='F').or.(UPLO=='f')) then
           call SGEMM('N','N',N,fpm(25),N,SONE,A,LDA,X(1,fpm(24)),N,SZERO,work(1,fpm(24)),N)
        else
           call SSYMM ('L', UPLO, N, fpm(25), SONE, A, LDA, X(1,fpm(24)), N, SZERO,work(1,fpm(24)), N)
        endif

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if ((UPLO=='F').or.(UPLO=='f')) then
           call SGEMM('N','N',N,fpm(25),N,SONE,B,LDB,X(1,fpm(24)),N,SZERO,work(1,fpm(24)),N)
        else
           call SSYMM ('L', UPLO, N, fpm(25), SONE, B, LDB, X(1,fpm(24)), N, SZERO,work(1,fpm(24)), N)
        endif

     end select
  end do



  call wdeallocate_2s(Aq)
  call wdeallocate_2s(Sq)
  call wdeallocate_2s(work)
  call wdeallocate_2c(workc)
  call wdeallocate_2c(Az)
  call wdeallocate_1i(ipivloc)


end subroutine sfeast_sygv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine sfeast_syev(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A REAL SYMMETRIC :: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) REAL SINGLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                        if option fpm(6)=1 selected                           
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2009-2012
  ! ====================================================================
  implicit none
  include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,LDA
  real,dimension(LDA,*):: A
  integer,dimension(*) :: fpm
  real :: epsout 
  integer :: loop
  real :: Emin,Emax
  integer :: M0
  real,dimension(*)  :: E
  real,dimension(N,*):: X
  integer :: mode
  real,dimension(*)    :: res
  integer :: info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc,i,s
  complex  :: Ze
  complex , dimension(:,:),pointer ::workc,Az
  real, dimension(:,:),pointer ::work,Aq,Sq
  integer, dimension(:),pointer ::ipivloc
  real,parameter :: SONE=1.0E0,SZERO=0.0E0
  complex ,parameter :: ONEC=(SONE,SZERO)

  INFO = 0
  IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
     INFO=-101
  ELSE IF( N<=0 ) THEN
     INFO = -102
  ELSE IF(LDA<N ) THEN
     INFO = -104
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'SFEAST_SYEV', -INFO+100 )
     RETURN
  END IF
  infoloc=0
  call wallocate_2s(Aq,M0,M0,infoloc)
  call wallocate_2s(Sq,M0,M0,infoloc)
  call wallocate_2s(work,N,M0,infoloc)
  call wallocate_2c(workc,N,M0,infoloc)
  call wallocate_2c(Az,N,N,infoloc)
  call wallocate_1i(ipivloc,N,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if


  ijob=-1 ! initialization
  do while (ijob/=0)
     call sfeast_srci(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)

     select case(ijob)
     case(10) !! factorize (zeB-A)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION

        if ((UPLO=='L').or.(UPLO=='l')) then
           do i=1,N-1
              s=N-i
              Az(i:N,i)=-A(i:N,i)*ONEC
              Az(i,i)=Az(i,i)+Ze
              call CCOPY(s,Az(i+1,i),1,Az(i,i+1),N)
           enddo
           Az(N,N)=Ze-A(N,N)*ONEC

        elseif ((UPLO=='U').or.(UPLO=='u')) then
           do i=1,N-1
              s=N-i
              Az(i,i:N)=-A(i,i:N)*ONEC
              Az(i,i)=Az(i,i)+Ze
              call CCOPY(s,Az(i,i+1),N,Az(i+1,i),1)
           enddo
           Az(N,N)=Ze-A(N,N)*ONEC


        else

           Az(1:N,1:N)=-A(1:N,1:N)*ONEC
           do i=1,N
              Az(i,i)=Az(i,i)+Ze
           enddo

        end if


        call CGETRF(N,N,Az,N,IPIVloc,INFOloc)     
        if (infoloc/=0) then
           info=-2
           return
        end if


     case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:M0) result in to workc

        call CGETRS( 'N', N, M0, Az, N, IPIVloc, workc, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if


     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if ((UPLO=='F').or.(UPLO=='f')) then
           call SGEMM('N','N',N,fpm(25),N,SONE,A,LDA,X(1,fpm(24)),N,SZERO,work(1,fpm(24)),N)
        else
           call SSYMM ('L', UPLO, N, fpm(25), SONE, A, LDA, X(1,fpm(24)), N, SZERO,work(1,fpm(24)), N)
        endif

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        call SLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )

     end select
  end do


  call wdeallocate_2s(Aq)
  call wdeallocate_2s(Sq)
  call wdeallocate_2s(work)
  call wdeallocate_2c(workc)
  call wdeallocate_2c(Az)
  call wdeallocate_1i(ipivloc)


end subroutine sfeast_syev



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine cfeast_hegv(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the generalized Ax=eBx eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
  !  B          (input)        COMPLEX SINGLE PRECISION (LDB,N):  Matrix B 
  !  LDB        (input)        INTEGER: Leading dimension of matrix B (LDB>=N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) COMPLEX SINGLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                        if option fpm(6)=1 selected                           
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2009-2012
  ! ====================================================================
  implicit none
  include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,LDA,LDB
  complex ,dimension(LDA,*):: A
  complex ,dimension(LDB,*):: B
  integer,dimension(*) :: fpm
  real :: epsout 
  integer :: loop
  real :: Emin,Emax
  integer :: M0
  real,dimension(*)  :: E
  complex ,dimension(N,*):: X
  integer :: mode
  real,dimension(*)    :: res
  integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc,i,s
  complex  :: Ze
  complex , dimension(:,:),pointer ::work1,work2,Az,zAq,zSq
  complex , dimension(:),pointer ::ztmp
  integer, dimension(:),pointer ::ipivloc
  real,parameter :: SONE=1.0E0,SZERO=0.0E0
  complex ,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)



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
     CALL XERBLA( 'CFEAST_HEGV', -INFO+100 )
     RETURN
  END IF


  infoloc=0
  call wallocate_2c(zAq,M0,M0,infoloc)
  call wallocate_2c(zSq,M0,M0,infoloc)
  call wallocate_2c(work1,N,M0,infoloc)
  call wallocate_2c(work2,N,M0,infoloc)
  call wallocate_2c(Az,N,N,infoloc)
  call wallocate_1i(ipivloc,N,infoloc)
  call wallocate_1c(ztmp,N,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if

  fpm(11)=1 ! half-contour only- requires solving (zB-A)^H x=f
  !fpm(11)=2 ! full contour  (2 half-contour)  - does not require solving for the transpose-conjg


  ijob=-1 ! initialization


  do while (ijob/=0) 
     call cfeast_hrci(ijob,N,Ze,work1,work2,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
     select case(ijob)

     case(10) !! factorize (zeB-A)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION

        if ((UPLO=='L').or.(UPLO=='l')) then
           do i=1,N-1
              s=N-i
              Az(i:N,i)=Ze*B(i:N,i)
              call CCOPY(s,Az(i+1,i),1,Az(i,i+1),N)
              call CLACGV(s, Az(i,i+1), N )
              call CSCAL(s,Ze/conjg(Ze),Az(i,i+1),N)
              Az(i:N,i)=Az(i:N,i)-A(i:N,i)
           enddo
           Az(N,N)=Ze*B(N,N)-A(N,N)*ONEC
           do i=1,N-1
              s=N-i
              call CCOPY(s,A(i+1,i),1,ztmp(1),1)
              call CLACGV(s,ztmp(1),1)
              call CAXPY(s,-ONEC,ztmp(1),1,Az(i,i+1),N)
           enddo


        elseif ((UPLO=='U').or.(UPLO=='u')) then

           do i=1,N-1
              s=N-i
              Az(i,i:N)=Ze*B(i,i:N)
              call CCOPY(s,Az(i,i+1),N,Az(i+1,i),1)
              call CLACGV(s, Az(i+1,i), 1 )
              call CSCAL(s,Ze/conjg(Ze),Az(i+1,i),1)
              Az(i,i:N)=Az(i,i:N)-A(i,i:N)
           enddo
           Az(N,N)=Ze*B(N,N)-A(N,N)*ONEC
           do i=1,N-1
              s=N-i
              call CCOPY(s,A(i,i+1),N,ztmp(1),1)
              call CLACGV(s,ztmp(1),1)
              call CAXPY(s,-ONEC,ztmp(1),1,Az(i+1,i),1)
           enddo

        else ! full 
           Az(1:N,1:N)=Ze*B(1:N,1:N)-A(1:N,1:N)
        end if




        call CGETRF(N,N,Az,N,IPIVloc,INFOloc)     
        if (infoloc/=0) then
           info=-2
           return
        end if

     case(11) !!solve the linear system (ZeB-A)x=work2(1:N,1:M0) result in to work2

        call CGETRS( 'N', N, M0, Az, N, IPIVloc, work2, N, INFOloc )

     case(21) !!solve the linear system (ZeB-A)^H x=work2(1:N,1:M0) result in to work2

        call CGETRS( 'C', N, M0, Az, N, IPIVloc, work2, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if

     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if ((UPLO=='F').or.(UPLO=='f')) then
           call CGEMM('N','N',N,fpm(25),N,ONEC,A,LDA,X(1,fpm(24)),N,ZEROC,work1(1,fpm(24)),N)
        else
           call CHEMM('L', UPLO, N, fpm(25), ONEC, A, LDA, X(1,fpm(24)), N, ZEROC,work1(1,fpm(24)), N)
        endif

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if ((UPLO=='F').or.(UPLO=='f')) then
           call CGEMM('N','N',N,fpm(25),N,ONEC,B,LDB,X(1,fpm(24)),N,ZEROC,work1(1,fpm(24)),N)
        else
           call CHEMM ('L', UPLO, N, fpm(25), ONEC, B, LDB, X(1,fpm(24)), N, ZEROC,work1(1,fpm(24)), N)
        endif

     end select
  end do

  call wdeallocate_2c(zAq)
  call wdeallocate_2c(zSq)
  call wdeallocate_2c(work1)
  call wdeallocate_2c(work2)
  call wdeallocate_2c(Az)
  call wdeallocate_1i(ipivloc)
  call wdeallocate_1c(ztmp)


end subroutine cfeast_hegv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine cfeast_heev(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
  !  Purpose 
  !  =======
  !  FEAST DENSE INTERFACE
  !  Solve the standard Ax=ex eigenvalue problem
  !  
  !  A COMPLEX HERMITIAN :: DENSE FORMAT 
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  UPLO       (input)       character: specifies whether the full part, or the upper or lower
  !                           triangular part of the matrix(ces) is(are) being supplied.
  !  N          (input)        INTEGER: Size system
  !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A 
  !  LDA        (input)        INTEGER: Leading dimension of matrix A (LDA>=N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! LIST of FEAST ARGUMENTS COMMON TO ALL FEAST INTERFACES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  fpm (input/output) INTEGER(*) : FEAST parameters
  !  epsout     (output)       REAL SINGLE PRECISION : Error on the trace
  !  loop       (output)       INTEGER : # of iterative loop to reach convergence 
  !  Emin,Emax  (input)        REAL SINGLE PRECISION: search interval
  !  M0         (input/output) INTEGER: Size subspace
  !  lambda     (output)       REAL SINGLE PRECISION(M0)   : Eigenvalues -solution
  !  q          (input/output) COMPLEX SINGLE PRECISION(N,M0) : 
  !                                                       On entry: subspace initial guess if fpm(5)=1 
  !                                                       On exit : Eigenvectors-solution
  !  mode       (output)       INTEGER : # of eigenvalues found in the search interval
  !  res        (output)       REAL SINGLE PRECISION(M0) : Relative Residual of the solution (1-norm)
  !                                                        if option fpm(6)=1 selected                           
  !  info       (output)       INTEGER: Error handling (0: successful exit)
  !=====================================================================
  ! Eric Polizzi 2009-2012
  ! ====================================================================
  implicit none
  include 'f90_noruntime_interface.fi'
  character(len=1) :: UPLO
  integer :: N,LDA
  complex ,dimension(LDA,*):: A
  integer,dimension(*) :: fpm
  real :: epsout  
  integer :: loop
  real :: Emin,Emax
  integer :: M0
  real,dimension(*)  :: E
  complex ,dimension(N,*):: X
  integer :: mode
  real,dimension(*)    :: res
  integer :: info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: ijob,infoloc,i,s
  complex  :: Ze
  complex , dimension(:,:),pointer ::work1,work2,Az,zAq,zSq
  integer, dimension(:),pointer ::ipivloc
  real,parameter :: SONE=1.0E0,SZERO=0.0E0
  complex ,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)


  INFO = 0
  IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
     INFO=-101
  ELSE IF( N<=0 ) THEN
     INFO = -102
  ELSE IF(LDA<N ) THEN
     INFO = -104
  END IF
  IF( INFO.NE.0 ) THEN
     CALL XERBLA( 'CFEAST_HEEV', -INFO+100 )
     RETURN
  END IF

  infoloc=0
  call wallocate_2c(zAq,M0,M0,infoloc)
  call wallocate_2c(zSq,M0,M0,infoloc)
  call wallocate_2c(work1,N,M0,infoloc)
  call wallocate_2c(work2,N,M0,infoloc)
  call wallocate_2c(Az,N,N,infoloc)
  call wallocate_1i(ipivloc,N,infoloc)
  if (infoloc/=0) then
     info=-1
     return
  end if

  fpm(11)=1 ! half-contour only- requires solving (zB-A)^H x=f
  !fpm(11)=2 ! full contour  (2 half-contour)  - does not require solving for the transpose-conjg


  ijob=-1 ! initialization


  do while (ijob/=0) 
     call cfeast_hrci(ijob,N,Ze,work1,work2,zAq,zSq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
     select case(ijob)
     case(10) !! factorize (zeB-A)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION

        if ((UPLO=='L').or.(UPLO=='l')) then
           do i=1,N-1
              s=N-i
              Az(i:N,i)=-A(i:N,i)
              Az(i,i)=Az(i,i)+Ze
              call CCOPY(s,Az(i+1,i),1,Az(i,i+1),N)
              call CLACGV(s, Az(i,i+1), N )
           enddo
           Az(N,N)=Ze-A(N,N)*ONEC

        elseif ((UPLO=='U').or.(UPLO=='u')) then
           do i=1,N-1
              s=N-i
              Az(i,i:N)=-A(i,i:N)
              Az(i,i)=Az(i,i)+Ze
              call CCOPY(s,Az(i,i+1),N,Az(i+1,i),1)
              call CLACGV(s, Az(i+1,i), 1 )
           enddo
           Az(N,N)=Ze-A(N,N)*ONEC

        else

           Az(1:N,1:N)=-A(1:N,1:N)
           do i=1,N
              Az(i,i)=Az(i,i)+Ze
           enddo


        end if


        call CGETRF(N,N,Az,N,IPIVloc,INFOloc)     
        if (infoloc/=0) then
           info=-2
           return
        end if


     case(11) !!solve the linear system (ZeB-A)x=work2(1:N,1:M0) result in to work2

        call CGETRS( 'N', N, M0, Az, N, IPIVloc, work2, N, INFOloc )
        if (infoloc/=0) then
           info=-2
           return
        end if


     case(21) !!solve the linear system (ZeB-A)^H x=work2(1:N,1:M0) result in to work2

        call CGETRS( 'C', N, M0, Az, N, IPIVloc, work2, N, INFOloc )


     case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        if ((UPLO=='F').or.(UPLO=='f')) then
           call CGEMM('N','N',N,fpm(25),N,ONEC,A,LDA,X(1,fpm(24)),N,ZEROC,work1(1,fpm(24)),N)
        else
           call CHEMM ('L', UPLO, N, fpm(25), ONEC, A, LDA, X(1,fpm(24)), N, ZEROC,work1(1,fpm(24)), N)
        endif

     case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

        call CLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work1(1,fpm(24)), N )

     end select
  end do

  call wdeallocate_2c(zAq)
  call wdeallocate_2c(zSq)
  call wdeallocate_2c(work1)
  call wdeallocate_2c(work2)
  call wdeallocate_2c(Az)
  call wdeallocate_1i(ipivloc)


end subroutine cfeast_heev


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






