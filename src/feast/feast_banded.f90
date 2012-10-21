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
!!!!!!!!!!!!!!!!! FEAST PREDEFINED BANDED INTERFACES !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! List of routines:
!-------------------

!{S,D,C,Z}FEAST_{SB,HB}{EV,GV}

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Symmetric eigenvalue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! single/double precision
! real symmetric 
!{s,d}feast_sbgv ! generalized
!{s,d}feast_sbev ! standard

! complex Hermitian
!{c,z}feast_hbgv ! generalized
!{c,z}feast_hbev ! standard
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  
  include 'lbprim.f90' !! Banded primitives
 
  subroutine dfeast_sbgv(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A
    !  LDA        (input)        INTEGER: Leading dimension of matrix A 
    !                                     If UPLO='F'  LDA>=2*kla+1
    !                                     If UPLO/='F' LDA>=kla+1
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  B          (input)        REAL DOUBLE PRECISION (LDB,N):  Matrix B
    !  LDB        (input)        INTEGER: Leading dimension of matrix B
    !                                     If UPLO='F'  LDB>=2*klb+1
    !                                     If UPLO/='F' LDB>=klb+1
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
    integer :: N,LDA,LDB,kla,klb
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
    integer :: ijob,i,band,klz,s,infoloc
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::workc,Az
    complex(kind=(kind(1.0d0))), dimension(:),pointer ::ztmp
    double precision, dimension(:,:),pointer ::work,Aq,Sq
    double precision :: norm,nzero
    double precision,parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)
    integer :: mlda,mldb
    mlda=2*kla+1
    mldb=2*klb+1
    IF ((UPLO=='L').or.(UPLO=='l').or.(UPLO=='U').or.(UPLO=='u')) THEN
       mlda=mlda-kla
       mldb=mldb-klb
    ENDIF

    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    ELSE IF (kla>=N) THEN
       INFO=-103
    ELSE IF(LDA<mlda) THEN
       INFO = -105
    ELSE IF (klb>=N) THEN
       INFO=-106
    ELSE IF(LDB<mldb) THEN
       INFO = -108
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DFEAST_SBGV', -INFO+100 )
       RETURN
    END IF


    klz=max(kla,klb)
    band=2*klz+1

    infoloc=0
    call wallocate_2d(Aq,M0,M0,infoloc)
    call wallocate_2d(Sq,M0,M0,infoloc)
    call wallocate_2d(work,N,M0,infoloc)
    call wallocate_2z(workc,N,M0,infoloc)
    call wallocate_2z(Az,band,N,infoloc)
    call wallocate_1z(ztmp,band,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if


    ijob=-1 ! initialization

    do while (ijob/=0)

       call dfeast_srci(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)

       select case(ijob)

       case(10) !! factorize (zeB-A)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION
          if (kla>=klb) then 

             if ((UPLO=='L').or.(UPLO=='l')) then
                Az(kla+1:band,1:N)=-A(1:kla+1,1:N)*ONEC
                Az(kla+1:kla+1+klb,1:N)=Az(kla+1:kla+1+klb,1:N)+Ze*B(1:klb+1,1:N)
                do i=1,N-1
                   s=min(klz,N-i)
                   call ZCOPY(s,Az(klz+2,i),1,Az(klz,i+1),2*klz)
                enddo

             elseif ((UPLO=='U').or.(UPLO=='u')) then
                Az(1:kla+1,1:N)=-A(1:kla+1,1:N)*ONEC
                Az(kla+1-klb:kla+1,1:N)=Az(kla+1-klb:kla+1,1:N)+Ze*B(1:klb+1,1:N)
                do i=1,N-1
                   s=min(klz,N-i)
                   call ZCOPY(s,Az(klz,i+1),2*klz,Az(klz+2,i),1)
                enddo

             else !UPLO=F

                do i=1,N
                   Az(1:band,i)=-A(1:band,i)*ONEC
                   ztmp(1:2*klb+1)=Ze*B(1:2*klb+1,i)
                   call ZAXPY(2*klb+1,ONEC,ztmp(1),1,Az(kla+1-klb,i),1)
                end do

             end if


          else ! kla<klb!!!!!!!!!!!!!!!!


             if ((UPLO=='L').or.(UPLO=='l')) then
                Az(klb+1:band,1:N)=Ze*B(1:klb+1,1:N)
                Az(klb+1:klb+1+kla,1:N)=Az(klb+1:klb+1+kla,1:N)-ONEC*A(1:kla+1,1:N)
                do i=1,N-1
                   s=min(klz,N-i)
                   call ZCOPY(s,Az(klz+2,i),1,Az(klz,i+1),2*klz)
                enddo

             elseif  ((UPLO=='U').or.(UPLO=='u')) then
                Az(1:klb+1,1:N)=Ze*B(1:klb+1,1:N)
                Az(klb+1-kla:klb+1,1:N)=Az(klb+1-kla:klb+1,1:N)-ONEC*A(1:kla+1,1:N)
                do i=1,N-1
                   s=min(klz,N-i)
                   call ZCOPY(s,Az(klz,i+1),2*klz,Az(klz+2,i),1)
                enddo

             else !UPLO=F

                do i=1,N
                   Az(1:band,i)=B(1:band,i)*Ze
                   ztmp(1:2*kla+1)=A(1:2*kla+1,i)
                   call ZAXPY(2*kla+1,-ONEC,ztmp(1),1,Az(klb+1-kla,i),1)
                end do


             endif

          end if




          nzero=DZERO
          call ZGBALU(N,klz,klz,Az,band,nzero,norm,infoloc)
          if (infoloc/=0) then
             info=-2
             return
          end if


       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:M0) result in to workc

          call ZTBSM('L','N','U',N,M0,klz,Az(klz+1,1),band,workc,N)
          call ZTBSM('U','N','N',N,M0,klz,Az(1,1),band,workc,N)


       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call DSBMM(UPLO,n,fpm(25),kla,DONE,A(1,1),LDA,X(1,fpm(24)),N,DZERO,work(1,fpm(24)),N)

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call DSBMM(UPLO,n,fpm(25),klb,DONE,B(1,1),LDB,X(1,fpm(24)),N,DZERO,work(1,fpm(24)),N)

       end select
    end do




    call wdeallocate_2d(Aq)
    call wdeallocate_2d(Sq)
    call wdeallocate_2z(workc)
    call wdeallocate_2z(Az)
    call wdeallocate_2d(work)
    call wdeallocate_1z(ztmp)

  end subroutine dfeast_sbgv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





  subroutine dfeast_sbev(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A REAL SYMMETRIC :: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        REAL DOUBLE PRECISION (LDA,N):  Matrix A
    !  LDA        (input)        INTEGER: Leading dimension of matrix A 
    !                                     If UPLO='F'  LDA>=2*kla+1
    !                                     If UPLO/='F' LDA>=kla+1
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
    integer :: N,LDA,kla
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
    integer :: ijob,i,band,klz,s,infoloc
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::workc,Az
    double precision, dimension(:,:),pointer ::work,Aq,Sq
    double precision :: norm,nzero
    double precision,parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO)
    integer :: mlda
    mlda=2*kla+1
    IF ((UPLO=='L').or.(UPLO=='l').or.(UPLO=='U').or.(UPLO=='u')) THEN
       mlda=mlda-kla
    ENDIF

    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    ELSE IF (kla>=N) THEN
       INFO=-103
    ELSE IF(LDA<mlda) THEN
       INFO = -105
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DFEAST_SBEV', -INFO+100 )
       RETURN
    END IF





    klz=kla
    band=2*klz+1

    infoloc=0
    call wallocate_2d(Aq,M0,M0,infoloc)
    call wallocate_2d(Sq,M0,M0,infoloc)
    call wallocate_2d(work,N,M0,infoloc)
    call wallocate_2z(workc,N,M0,infoloc)
    call wallocate_2z(Az,band,N,infoloc)



    ijob=-1 ! initialization

    do while (ijob/=0)

       call dfeast_srci(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)


       select case(ijob)

       case(10) !! factorize (zeB-A)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION
          if ((UPLO=='L').or.(UPLO=='l')) then
             Az(klz+1:band,1:N)=-A(1:kla+1,1:N)*ONEC
             do i=1,N-1
                s=min(klz,N-i) 
                call ZCOPY(s,Az(klz+2,i),1,Az(klz,i+1),2*klz)
                Az(klz+1,i)=-A(1,i)+Ze
             enddo
             Az(klz+1,N)=-A(1,N)+Ze


          elseif ((UPLO=='U').or.(UPLO=='u')) then
             Az(1:klz+1,1:N)=-A(1:klz+1,1:N)*ONEC
             do i=1,N-1
                s=min(klz,N-i)
                call ZCOPY(s,Az(klz,i+1),2*klz,Az(klz+2,i),1)
                Az(klz+1,i)=-A(klz+1,i)+Ze
             enddo
             Az(klz+1,N)=-A(klz+1,N)+Ze


          else !UPLO=F

             do i=1,N
                Az(1:band,i)=-A(1:band,i)*ONEC
                Az(klz+1,i)=Az(klz+1,i)+Ze
             end do

          end if

          nzero=DZERO
          call ZGBALU(N,klz,klz,Az,band,nzero,norm,infoloc)
          if (infoloc/=0) then
             info=-2
             return
          end if


       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:M0) result in to workc

          call ZTBSM('L','N','U',N,M0,klz,Az(klz+1,1),band,workc,N)
          call ZTBSM('U','N','N',N,M0,klz,Az(1,1),band,workc,N)


       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call DSBMM(UPLO,n,fpm(25),kla,DONE,A(1,1),LDA,X(1,fpm(24)),N,DZERO,work(1,fpm(24)),N)

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call DLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )

       end select
    end do




    call wdeallocate_2d(Aq)
    call wdeallocate_2d(Sq)
    call wdeallocate_2z(workc)
    call wdeallocate_2z(Az)
    call wdeallocate_2d(work)

  end subroutine dfeast_sbev


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine zfeast_hbgv(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A
    !  LDA        (input)        INTEGER: Leading dimension of matrix A 
    !                                     If UPLO='F'  LDA>=2*kla+1
    !                                     If UPLO/='F' LDA>=kla+1
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  B          (input)        COMPLEX DOUBLE PRECISION (LDB,N):  Matrix B
    !  LDB        (input)        INTEGER: Leading dimension of matrix B
    !                                     If UPLO='F'  LDB>=2*klb+1
    !                                     If UPLO/='F' LDB>=klb+1
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
    integer :: N,LDA,LDB,kla,klb
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,i,band,klz,s,infoloc
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::work1,work2,Az,zAq,zSq
    complex(kind=(kind(1.0d0))), dimension(:),pointer ::ztmp
    double precision :: norm,nzero
    double precision,parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
    integer :: mlda,mldb
    mlda=2*kla+1
    mldb=2*klb+1
    IF ((UPLO=='L').or.(UPLO=='l').or.(UPLO=='U').or.(UPLO=='u')) THEN
       mlda=mlda-kla
       mldb=mldb-klb
    ENDIF



    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    ELSE IF (kla>=N) THEN
       INFO=-103
    ELSE IF(LDA<mlda) THEN
       INFO = -105
    ELSE IF (klb>=N) THEN
       INFO=-106
    ELSE IF(LDB<mldb) THEN
       INFO = -108
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'ZFEAST_HBGV', -INFO+100 )
       RETURN
    END IF



    klz=max(kla,klb)
    band=2*klz+1

    infoloc=0
    call wallocate_2z(zAq,M0,M0,infoloc)
    call wallocate_2z(zSq,M0,M0,infoloc)
    call wallocate_2z(work1,N,M0,infoloc)
    call wallocate_2z(work2,N,M0,infoloc)
    call wallocate_2z(Az,band,N,infoloc)
    call wallocate_1z(ztmp,klz,infoloc)
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
          if (kla>=klb) then 

             if ((UPLO=='L').or.(UPLO=='l')) then

                Az(kla+1:band,1:N)=-A(1:kla+1,1:N)
                do i=1,N-1
                   s=min(klz,N-i)
                   call ZCOPY(s,Az(klz+2,i),1,Az(klz,i+1),2*klz)
                   call ZLACGV(s, Az(klz,i+1), 2*klz )
                enddo
                Az(kla+1:kla+1+klb,1:N)=Az(kla+1:kla+1+klb,1:N)+Ze*B(1:klb+1,1:N)
                do i=1,N-1
                   s=min(klb,N-i)
                   call ZCOPY(s,B(2,i),1,ztmp(1),1)
                   call ZLACGV(s,ztmp(1),1)
                   call ZAXPY(s,Ze,ztmp(1),1,Az(klz,i+1),2*klz)
                enddo


             elseif ((UPLO=='U').or.(UPLO=='u')) then


                Az(1:kla+1,1:N)=-A(1:kla+1,1:N)
                do i=1,N-1
                   s=min(klz,N-i)
                   call ZCOPY(s,Az(klz,i+1),2*klz,Az(klz+2,i),1)
                   call ZLACGV(s, Az(klz+2,i), 1 )
                enddo
                Az(kla+1-klb:kla+1,1:N)=Az(kla+1-klb:kla+1,1:N)+Ze*B(1:klb+1,1:N)
                do i=1,N-1
                   s=min(klb,N-i)
                   call ZCOPY(s,B(klb,i+1),LDB-1,ztmp(1),1)
                   call ZLACGV(s,ztmp(1),1)
                   call ZAXPY(s,Ze,ztmp(1),1,Az(klz+2,i),1)
                enddo



             else !UPLO=F


                do i=1,N
                   call ZCOPY(band,A(1,i),1,Az(1,i),1)
                   call ZSCAL(band,-ONEC,Az(1,i),1)
                   call ZAXPY(2*klb+1,Ze,B(1,i),1,Az(kla+1-klb,i),1)
                end do

             end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          else ! kla<klb!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


             if ((UPLO=='L').or.(UPLO=='l')) then


                Az(klz+1:band,1:N)=Ze*B(1:klb+1,1:N)
                do i=1,N-1
                   s=min(klz,N-i)
                   call ZCOPY(s,Az(klz+2,i),1,Az(klz,i+1),2*klz)
                   call ZLACGV(s, Az(klz,i+1), 2*klz )
                   call ZSCAL(s,Ze/conjg(Ze),Az(klz,i+1),2*klz)
                enddo
                Az(klz+1:klz+1+kla,1:N)=Az(klb+1:klb+1+kla,1:N)-A(1:kla+1,1:N)
                do i=1,N-1
                   s=min(kla,N-i)
                   call ZCOPY(s,A(2,i),1,ztmp(1),1)
                   call ZLACGV(s,ztmp(1),1)
                   call ZAXPY(s,-ONEC,ztmp(1),1,Az(klz,i+1),2*klz)
                enddo

             elseif  ((UPLO=='U').or.(UPLO=='u')) then


                Az(1:klz+1,1:N)=Ze*B(1:klb+1,1:N)
                do i=1,N-1
                   s=min(klz,N-i)
                   call ZCOPY(s,Az(klz,i+1),2*klz,Az(klz+2,i),1)
                   call ZLACGV(s, Az(klz+2,i), 1 )
                   call ZSCAL(s,Ze/conjg(Ze),Az(klz+2,i),1)
                enddo
                Az(klz+1-kla:klz+1,1:N)=Az(klz+1-kla:klz+1,1:N)-A(1:kla+1,1:N)
                do i=1,N-1
                   s=min(kla,N-i)
                   call ZCOPY(s,A(kla,i+1),LDA-1,ztmp(1),1)
                   call ZLACGV(s,ztmp(1),1)
                   call ZAXPY(s,-ONEC,ztmp(1),1,Az(klz+2,i),1)
                enddo



             else !UPLO=F


                do i=1,N
                   call ZCOPY(band,B(1,i),1,Az(1,i),1)
                   call ZSCAL(band,Ze,Az(1,i),1)
                   call ZAXPY(2*kla+1,-ONEC,A(1,i),1,Az(klb+1-kla,i),1)
                end do


             endif

          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






          nzero=DZERO
          call ZGBALU(N,klz,klz,Az,band,nzero,norm,infoloc)
          if (infoloc/=0) then
             info=-2
             return
          end if


       case(11) !!solve the linear system (ZeB-A)x=work2(1:N,1:M0) result in to work2

          call ZTBSM('L','N','U',N,M0,klz,Az(klz+1,1),band,work2,N)
          call ZTBSM('U','N','N',N,M0,klz,Az(1,1),band,work2,N)


       case(21) !!solve the linear system (ZeB-A)^H x=work2(1:N,1:M0) result in to work2

          call ZTBSM('U','C','N',N,M0,klz,Az(1,1),band,work2,N)
          call ZTBSM('L','C','U',N,M0,klz,Az(klz+1,1),band,work2,N)



       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call ZHBMM(UPLO,n,fpm(25),kla,ONEC,A(1,1),LDA,X(1,fpm(24)),N,ZEROC,work1(1,fpm(24)),N)

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call ZHBMM(UPLO,n,fpm(25),klb,ONEC,B(1,1),LDB,X(1,fpm(24)),N,ZEROC,work1(1,fpm(24)),N)


       end select
    end do




    call wdeallocate_2z(zAq)
    call wdeallocate_2z(zSq)
    call wdeallocate_2z(work2)
    call wdeallocate_2z(Az)
    call wdeallocate_2z(work1)
    call wdeallocate_1z(ztmp)


  end subroutine zfeast_hbgv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine zfeast_hbev(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: BANDED FORMAT 
    ! 
    !  DOUBLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX DOUBLE PRECISION (LDA,N):  Matrix A
    !  LDA        (input)        INTEGER: Leading dimension of matrix A 
    !                                     If UPLO='F'  LDA>=2*kla+1
    !                                     If UPLO/='F' LDA>=kla+1
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
    integer :: N,LDA,kla
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: ijob,i,band,klz,s,infoloc
    complex(kind=(kind(1.0d0))) :: Ze
    complex(kind=(kind(1.0d0))), dimension(:,:),pointer ::work1,work2,Az,zAq,zSq
    double precision :: norm,nzero
    double precision,parameter :: DONE=1.0d0, DZERO=0.0d0
    complex(kind=(kind(1.0d0))),parameter :: ONEC=(DONE,DZERO),ZEROC=(DZERO,DZERO)
    integer :: mlda
    mlda=2*kla+1
    IF ((UPLO=='L').or.(UPLO=='l').or.(UPLO=='U').or.(UPLO=='u')) THEN
       mlda=mlda-kla
    ENDIF


    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    ELSE IF (kla>=N) THEN
       INFO=-103
    ELSE IF(LDA<mlda) THEN
       INFO = -105
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'ZFEAST_HBEV', -INFO+100 )
       RETURN
    END IF



    klz=kla
    band=2*klz+1

    infoloc=0
    call wallocate_2z(zAq,M0,M0,infoloc)
    call wallocate_2z(zSq,M0,M0,infoloc)
    call wallocate_2z(work1,N,M0,infoloc)
    call wallocate_2z(work2,N,M0,infoloc)
    call wallocate_2z(Az,band,N,infoloc)
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

          !Az(1:band,1:N)=-A(1:band,1:N)
          !do i=1,N
          !Az(klz+1,i)=Az(klz+1,i)+Ze
          !enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION
          if ((UPLO=='L').or.(UPLO=='l')) then

             do i=1,N-1
                s=min(klz,N-i)
                call ZCOPY(s,A(2,i),1,Az(klz+2,i),1)
                call ZSCAL(s,-ONEC,Az(klz+2,i),1)
                call ZCOPY(s,Az(klz+2,i),1,Az(klz,i+1),2*klz)
                call ZLACGV(s, Az(klz,i+1), 2*klz )
                Az(klz+1,i)=-A(1,i)+Ze
             enddo
             Az(klz+1,N)=-A(1,N)+Ze

          elseif ((UPLO=='U').or.(UPLO=='u')) then

             do i=1,N-1
                s=min(klz,N-i)
                call ZCOPY(s,A(klz,i+1),LDA-1,Az(klz,i+1),2*klz)
                call ZSCAL(s,-ONEC,Az(klz,i+1),2*klz)
                call ZCOPY(s,Az(klz,i+1),2*klz,Az(klz+2,i),1)
                call ZLACGV(s, Az(klz+2,i),1)
                Az(klz+1,i)=-A(klz+1,i)+Ze
             enddo
             Az(klz+1,N)=-A(klz+1,N)+Ze


          else !UPLO=F

             do i=1,N
                call ZCOPY(band,A(1,i),1,Az(1,i),1)
                call ZSCAL(band,-ONEC,Az(1,i),1)
                Az(klz+1,i)=Az(klz+1,i)+Ze
             end do


          end if



          nzero=DZERO
          call ZGBALU(N,klz,klz,Az,band,nzero,norm,infoloc)
          if (infoloc/=0) then
             info=-2
             return
          end if


       case(11) !!solve the linear system (ZeB-A)x=work2(1:N,1:M0) result in to work2



          call ZTBSM('L','N','U',N,M0,klz,Az(klz+1,1),band,work2,N)
          call ZTBSM('U','N','N',N,M0,klz,Az(1,1),band,work2,N)




       case(21) !!solve the linear system (ZeB-A)^H x=work2(1:N,1:M0) result in to work2

          call ZTBSM('U','C','N',N,M0,klz,Az(1,1),band,work2,N)
          call ZTBSM('L','C','U',N,M0,klz,Az(klz+1,1),band,work2,N)


       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call ZHBMM(UPLO,n,fpm(25),kla,ONEC,A(1,1),LDA,X(1,fpm(24)),N,ZEROC,work1(1,fpm(24)),N)

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call ZLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work1(1,fpm(24)), N )


       end select
    end do



    call wdeallocate_2z(zAq)
    call wdeallocate_2z(zSq)
    call wdeallocate_2z(work2)
    call wdeallocate_2z(Az)
    call wdeallocate_2z(work1)

  end subroutine zfeast_hbev
 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine sfeast_sbgv(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A REAL SYMMETRIC, B SYMMETRIC POSITIVE DEFINITE:: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A
    !  LDA        (input)        INTEGER: Leading dimension of matrix A 
    !                                     If UPLO='F'  LDA>=2*kla+1
    !                                     If UPLO/='F' LDA>=kla+1
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  B          (input)        REAL SINGLE PRECISION (LDB,N):  Matrix B
    !  LDB        (input)        INTEGER: Leading dimension of matrix B
    !                                     If UPLO='F'  LDB>=2*klb+1
    !                                     If UPLO/='F' LDB>=klb+1
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
    integer :: N,LDA,LDB,kla,klb
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
    integer :: ijob,i,band,klz,s,infoloc
    complex  :: Ze
    complex , dimension(:,:),pointer ::workc,Az
    complex , dimension(:),pointer ::ztmp
    real, dimension(:,:),pointer ::work,Aq,Sq
    real :: norm,nzero
    real,parameter :: SONE=1.0E0, SZERO=0.0E0
    complex ,parameter :: ONEC=(SONE,SZERO)
    integer :: mlda,mldb
    mlda=2*kla+1
    mldb=2*klb+1
    IF ((UPLO=='L').or.(UPLO=='l').or.(UPLO=='U').or.(UPLO=='u')) THEN
       mlda=mlda-kla
       mldb=mldb-klb
    ENDIF



    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    ELSE IF (kla>=N) THEN
       INFO=-103
    ELSE IF(LDA<mlda) THEN
       INFO = -105
    ELSE IF (klb>=N) THEN
       INFO=-106
    ELSE IF(LDB<mldb ) THEN
       INFO = -108
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'SFEAST_SBGV', -INFO+100 )
       RETURN
    END IF


    klz=max(kla,klb)
    band=2*klz+1

    infoloc=0
    call wallocate_2s(Aq,M0,M0,infoloc)
    call wallocate_2s(Sq,M0,M0,infoloc)
    call wallocate_2s(work,N,M0,infoloc)
    call wallocate_2c(workc,N,M0,infoloc)
    call wallocate_2c(Az,band,N,infoloc)
    call wallocate_1c(ztmp,band,infoloc)
    if (infoloc/=0) then
       info=-1
       return
    end if


    ijob=-1 ! initialization

    do while (ijob/=0)

       call sfeast_srci(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)

       select case(ijob)

       case(10) !! factorize (zeB-A)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION
          if (kla>=klb) then 

             if ((UPLO=='L').or.(UPLO=='l')) then
                Az(kla+1:band,1:N)=-A(1:kla+1,1:N)*ONEC
                Az(kla+1:kla+1+klb,1:N)=Az(kla+1:kla+1+klb,1:N)+Ze*B(1:klb+1,1:N)
                do i=1,N-1
                   s=min(klz,N-i)
                   call CCOPY(s,Az(klz+2,i),1,Az(klz,i+1),2*klz)
                enddo

             elseif ((UPLO=='U').or.(UPLO=='u')) then
                Az(1:kla+1,1:N)=-A(1:kla+1,1:N)*ONEC
                Az(kla+1-klb:kla+1,1:N)=Az(kla+1-klb:kla+1,1:N)+Ze*B(1:klb+1,1:N)
                do i=1,N-1
                   s=min(klz,N-i)
                   call CCOPY(s,Az(klz,i+1),2*klz,Az(klz+2,i),1)
                enddo

             else !UPLO=F

                do i=1,N
                   Az(1:band,i)=-A(1:band,i)*ONEC
                   ztmp(1:2*klb+1)=Ze*B(1:2*klb+1,i)
                   call CAXPY(2*klb+1,ONEC,ztmp(1),1,Az(kla+1-klb,i),1)
                end do

             end if


          else ! kla<klb!!!!!!!!!!!!!!!!


             if ((UPLO=='L').or.(UPLO=='l')) then
                Az(klb+1:band,1:N)=Ze*B(1:klb+1,1:N)
                Az(klb+1:klb+1+kla,1:N)=Az(klb+1:klb+1+kla,1:N)-ONEC*A(1:kla+1,1:N)
                do i=1,N-1
                   s=min(klz,N-i)
                   call CCOPY(s,Az(klz+2,i),1,Az(klz,i+1),2*klz)
                enddo

             elseif  ((UPLO=='U').or.(UPLO=='u')) then
                Az(1:klb+1,1:N)=Ze*B(1:klb+1,1:N)
                Az(klb+1-kla:klb+1,1:N)=Az(klb+1-kla:klb+1,1:N)-ONEC*A(1:kla+1,1:N)
                do i=1,N-1
                   s=min(klz,N-i)
                   call CCOPY(s,Az(klz,i+1),2*klz,Az(klz+2,i),1)
                enddo

             else !UPLO=F


                do i=1,N
                   Az(1:band,i)=B(1:band,i)*Ze
                   ztmp(1:2*kla+1)=A(1:2*kla+1,i)
                   call CAXPY(2*kla+1,-ONEC,ztmp(1),1,Az(klb+1-kla,i),1)
                end do

             endif
          end if


          nzero=SZERO
          call CGBALU(N,klz,klz,Az,band,nzero,norm,infoloc)
          if (infoloc/=0) then
             info=-2
             return
          end if


       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:M0) result in to workc

          call CTBSM('L','N','U',N,M0,klz,Az(klz+1,1),band,workc,N)
          call CTBSM('U','N','N',N,M0,klz,Az(1,1),band,workc,N)


       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)


          call SSBMM(UPLO,n,fpm(25),kla,SONE,A(1,1),LDA,X(1,fpm(24)),N,SZERO,work(1,fpm(24)),N)

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call SSBMM(UPLO,n,fpm(25),klb,SONE,B(1,1),LDB,X(1,fpm(24)),N,SZERO,work(1,fpm(24)),N)

       end select
    end do




    call wdeallocate_2s(Aq)
    call wdeallocate_2s(Sq)
    call wdeallocate_2c(workc)
    call wdeallocate_2c(Az)
    call wdeallocate_2s(work)
    call wdeallocate_1c(ztmp)

  end subroutine sfeast_sbgv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





  subroutine sfeast_sbev(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A REAL SYMMETRIC :: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        REAL SINGLE PRECISION (LDA,N):  Matrix A
    !  LDA        (input)        INTEGER: Leading dimension of matrix A 
    !                                     If UPLO='F'  LDA>=2*kla+1
    !                                     If UPLO/='F' LDA>=kla+1
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
    integer :: N,LDA,kla
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
    integer :: ijob,i,band,klz,s,infoloc
    complex  :: Ze
    complex , dimension(:,:),pointer ::workc,Az
    real, dimension(:,:),pointer ::work,Aq,Sq
    real :: norm,nzero
    real,parameter :: SONE=1.0E0, SZERO=0.0E0
    complex ,parameter :: ONEC=(SONE,SZERO)
    integer :: mlda
    mlda=2*kla+1
    IF ((UPLO=='L').or.(UPLO=='l').or.(UPLO=='U').or.(UPLO=='u')) THEN
       mlda=mlda-kla
    ENDIF



    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    ELSE IF (kla>=N) THEN
       INFO=-103
    ELSE IF(LDA<mlda) THEN
       INFO = -105
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'SFEAST_SBEV', -INFO+100 )
       RETURN
    END IF


    klz=kla
    band=2*klz+1

    infoloc=0
    call wallocate_2s(Aq,M0,M0,infoloc)
    call wallocate_2s(Sq,M0,M0,infoloc)
    call wallocate_2s(work,N,M0,infoloc)
    call wallocate_2c(workc,N,M0,infoloc)
    call wallocate_2c(Az,band,N,infoloc)



    ijob=-1 ! initialization

    do while (ijob/=0)

       call sfeast_srci(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)


       select case(ijob)

       case(10) !! factorize (zeB-A)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Format CONVERSION
          if ((UPLO=='L').or.(UPLO=='l')) then
             Az(klz+1:band,1:N)=-A(1:kla+1,1:N)*ONEC
             do i=1,N-1
                s=min(klz,N-i) 
                call CCOPY(s,Az(klz+2,i),1,Az(klz,i+1),2*klz)
                Az(klz+1,i)=-A(1,i)+Ze
             enddo
             Az(klz+1,N)=-A(1,N)+Ze


          elseif ((UPLO=='U').or.(UPLO=='u')) then
             Az(1:klz+1,1:N)=-A(1:klz+1,1:N)*ONEC
             do i=1,N-1
                s=min(klz,N-i)
                call CCOPY(s,Az(klz,i+1),2*klz,Az(klz+2,i),1)
                Az(klz+1,i)=-A(klz+1,i)+Ze
             enddo
             Az(klz+1,N)=-A(klz+1,N)+Ze


          else !UPLO=F

             do i=1,N
                Az(1:band,i)=-A(1:band,i)*ONEC
                Az(klz+1,i)=Az(klz+1,i)+Ze
             end do

          end if


          nzero=SZERO
          call CGBALU(N,klz,klz,Az,band,nzero,norm,infoloc)
          if (infoloc/=0) then
             info=-2
             return
          end if


       case(11) !!solve the linear system (ZeB-A)x=workc(1:N,1:M0) result in to workc

          call CTBSM('L','N','U',N,M0,klz,Az(klz+1,1),band,workc,N)
          call CTBSM('U','N','N',N,M0,klz,Az(1,1),band,workc,N)


       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call SSBMM(UPLO,n,fpm(25),kla,SONE,A(1,1),LDA,X(1,fpm(24)),N,SZERO,work(1,fpm(24)),N)

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call SLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work(1,fpm(24)), N )

       end select
    end do




    call wdeallocate_2s(Aq)
    call wdeallocate_2s(Sq)
    call wdeallocate_2c(workc)
    call wdeallocate_2c(Az)
    call wdeallocate_2s(work)

  end subroutine sfeast_sbev


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cfeast_hbgv(UPLO,N,kla,A,LDA,klb,B,LDB,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the generalized Ax=eBx eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN, B HERMITIAN POSITIVE DEFINITE:: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A
    !  LDA        (input)        INTEGER: Leading dimension of matrix A 
    !                                     If UPLO='F'  LDA>=2*kla+1
    !                                     If UPLO/='F' LDA>=kla+1
    !  klb        (input)        INTEGER: # of subdiagonals within the band of B 
    !  B          (input)        COMPLEX SINGLE PRECISION (LDB,N):  Matrix B
    !  LDB        (input)        INTEGER: Leading dimension of matrix B
    !                                     If UPLO='F'  LDB>=2*klb+1
    !                                     If UPLO/='F' LDB>=klb+1
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
    integer :: N,LDA,LDB,kla,klb
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
    integer :: ijob,i,band,klz,s,infoloc
    complex  :: Ze
    complex , dimension(:,:),pointer ::work1,work2,Az,zAq,zSq
    complex , dimension(:),pointer ::ztmp
    real :: norm,nzero
    real,parameter :: SONE=1.0E0, SZERO=0.0E0
    complex ,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)
    integer :: mlda,mldb
    mlda=2*kla+1
    mldb=2*klb+1
    IF ((UPLO=='L').or.(UPLO=='l').or.(UPLO=='U').or.(UPLO=='u')) THEN
       mlda=mlda-kla
       mldb=mldb-klb
    ENDIF




    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    ELSE IF (kla>=N) THEN
       INFO=-103
    ELSE IF(LDA<mlda) THEN
       INFO = -105
    ELSE IF (klb>=N) THEN
       INFO=-106
    ELSE IF(LDB<mldb) THEN
       INFO = -108
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'CFEAST_HBGV', -INFO+100 )
       RETURN
    END IF

    klz=max(kla,klb)
    band=2*klz+1

    infoloc=0
    call wallocate_2c(zAq,M0,M0,infoloc)
    call wallocate_2c(zSq,M0,M0,infoloc)
    call wallocate_2c(work1,N,M0,infoloc)
    call wallocate_2c(work2,N,M0,infoloc)
    call wallocate_2c(Az,band,N,infoloc)
    call wallocate_1c(ztmp,klz,infoloc)
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
          if (kla>=klb) then 

             if ((UPLO=='L').or.(UPLO=='l')) then

                Az(kla+1:band,1:N)=-A(1:kla+1,1:N)
                do i=1,N-1
                   s=min(klz,N-i)
                   call CCOPY(s,Az(klz+2,i),1,Az(klz,i+1),2*klz)
                   call CLACGV(s, Az(klz,i+1), 2*klz )
                enddo
                Az(kla+1:kla+1+klb,1:N)=Az(kla+1:kla+1+klb,1:N)+Ze*B(1:klb+1,1:N)
                do i=1,N-1
                   s=min(klb,N-i)
                   call CCOPY(s,B(2,i),1,ztmp(1),1)
                   call CLACGV(s,ztmp(1),1)
                   call CAXPY(s,Ze,ztmp(1),1,Az(klz,i+1),2*klz)
                enddo


             elseif ((UPLO=='U').or.(UPLO=='u')) then

                Az(1:kla+1,1:N)=-A(1:kla+1,1:N)
                do i=1,N-1
                   s=min(klz,N-i)
                   call CCOPY(s,Az(klz,i+1),2*klz,Az(klz+2,i),1)
                   call CLACGV(s, Az(klz+2,i), 1 )
                enddo
                Az(kla+1-klb:kla+1,1:N)=Az(kla+1-klb:kla+1,1:N)+Ze*B(1:klb+1,1:N)
                do i=1,N-1
                   s=min(klb,N-i)
                   call CCOPY(s,B(klb,i+1),LDB-1,ztmp(1),1)
                   call CLACGV(s,ztmp(1),1)
                   call CAXPY(s,Ze,ztmp(1),1,Az(klz+2,i),1)
                enddo



             else !UPLO=F

                do i=1,N
                   call CCOPY(band,A(1,i),1,Az(1,i),1)
                   call CSCAL(band,-ONEC,Az(1,i),1)
                   call CAXPY(2*klb+1,Ze,B(1,i),1,Az(kla+1-klb,i),1)
                end do

             end if


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          else ! kla<klb!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


             if ((UPLO=='L').or.(UPLO=='l')) then

                Az(klz+1:band,1:N)=Ze*B(1:klb+1,1:N)
                do i=1,N-1
                   s=min(klz,N-i)
                   call CCOPY(s,Az(klz+2,i),1,Az(klz,i+1),2*klz)
                   call CLACGV(s, Az(klz,i+1), 2*klz )
                   call CSCAL(s,Ze/conjg(Ze),Az(klz,i+1),2*klz)
                enddo
                Az(klz+1:klz+1+kla,1:N)=Az(klb+1:klb+1+kla,1:N)-A(1:kla+1,1:N)
                do i=1,N-1
                   s=min(kla,N-i)
                   call CCOPY(s,A(2,i),1,ztmp(1),1)
                   call CLACGV(s,ztmp(1),1)
                   call CAXPY(s,-ONEC,ztmp(1),1,Az(klz,i+1),2*klz)
                enddo

             elseif  ((UPLO=='U').or.(UPLO=='u')) then


                Az(1:klz+1,1:N)=Ze*B(1:klb+1,1:N)
                do i=1,N-1
                   s=min(klz,N-i)
                   call CCOPY(s,Az(klz,i+1),2*klz,Az(klz+2,i),1)
                   call CLACGV(s, Az(klz+2,i), 1 )
                   call CSCAL(s,Ze/conjg(Ze),Az(klz+2,i),1)
                enddo
                Az(klz+1-kla:klz+1,1:N)=Az(klz+1-kla:klz+1,1:N)-A(1:kla+1,1:N)
                do i=1,N-1
                   s=min(kla,N-i)
                   call CCOPY(s,A(kla,i+1),LDA-1,ztmp(1),1)
                   call CLACGV(s,ztmp(1),1)
                   call CAXPY(s,-ONEC,ztmp(1),1,Az(klz+2,i),1)
                enddo



             else !UPLO=F


                do i=1,N
                   call CCOPY(band,B(1,i),1,Az(1,i),1)
                   call CSCAL(band,Ze,Az(1,i),1)
                   call CAXPY(2*kla+1,-ONEC,A(1,i),1,Az(klb+1-kla,i),1)
                end do


             endif

          end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






          nzero=SZERO
          call CGBALU(N,klz,klz,Az,band,nzero,norm,infoloc)
          if (infoloc/=0) then
             info=-2
             return
          end if


       case(11) !!solve the linear system (ZeB-A)x=work2(1:N,1:M0) result in to work2

          call CTBSM('L','N','U',N,M0,klz,Az(klz+1,1),band,work2,N)
          call CTBSM('U','N','N',N,M0,klz,Az(1,1),band,work2,N)


       case(21) !!solve the linear system (ZeB-A)^H x=work2(1:N,1:M0) result in to work2

          call CTBSM('U','C','N',N,M0,klz,Az(1,1),band,work2,N)
          call CTBSM('L','C','U',N,M0,klz,Az(klz+1,1),band,work2,N)


       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call CHBMM(UPLO,n,fpm(25),kla,ONEC,A(1,1),LDA,X(1,fpm(24)),N,ZEROC,work1(1,fpm(24)),N)

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call CHBMM(UPLO,n,fpm(25),klb,ONEC,B(1,1),LDB,X(1,fpm(24)),N,ZEROC,work1(1,fpm(24)),N)


       end select
    end do




    call wdeallocate_2c(zAq)
    call wdeallocate_2c(zSq)
    call wdeallocate_2c(work2)
    call wdeallocate_2c(Az)
    call wdeallocate_2c(work1)
    call wdeallocate_1c(ztmp)


  end subroutine cfeast_hbgv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine cfeast_hbev(UPLO,N,kla,A,LDA,fpm,epsout,loop,Emin,Emax,M0,E,X,mode,res,info)
    !  Purpose 
    !  =======
    !  FEAST BANDED INTERFACE
    !  Solve the standard Ax=ex eigenvalue problem
    !  
    !  A COMPLEX HERMITIAN :: BANDED FORMAT 
    ! 
    !  SINGLE PRECISION version  
    !
    !  Arguments
    !  =========
    !
    !  UPLO       (input)        CHARACTER: specifies whether the full part, or the upper or lower
    !                                       triangular part of the matrix(ces) is(are) being supplied.
    !  N          (input)        INTEGER: Size system
    !  kla        (input)        INTEGER: # of subdiagonals within the band of A 
    !  A          (input)        COMPLEX SINGLE PRECISION (LDA,N):  Matrix A
    !  LDA        (input)        INTEGER: Leading dimension of matrix A 
    !                                     If UPLO='F'  LDA>=2*kla+1
    !                                     If UPLO/='F' LDA>=kla+1
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
    integer :: N,LDA,kla
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
    integer :: ijob,i,band,klz,s,infoloc
    complex  :: Ze
    complex , dimension(:,:),pointer ::work1,work2,Az,zAq,zSq
    real :: norm,nzero
    real,parameter :: SONE=1.0E0, SZERO=0.0E0
    complex ,parameter :: ONEC=(SONE,SZERO),ZEROC=(SZERO,SZERO)
    integer :: mlda
    mlda=2*kla+1
    IF ((UPLO=='L').or.(UPLO=='l').or.(UPLO=='U').or.(UPLO=='u')) THEN
       mlda=mlda-kla
    ENDIF



    INFO = 0
    IF ((UPLO/='F').and.(UPLO/='f').and.(UPLO/='L').and.(UPLO/='l').and.(UPLO/='U').and.(UPLO/='u')) THEN
       INFO=-101
    ELSE IF ( N<=0 ) THEN
       INFO = -102
    ELSE IF (kla>=N) THEN
       INFO=-103
    ELSE IF(LDA<mlda) THEN
       INFO = -105
    END IF
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'CFEAST_HBEV', -INFO+100 )
       RETURN
    END IF



    klz=kla
    band=2*klz+1

    infoloc=0
    call wallocate_2c(zAq,M0,M0,infoloc)
    call wallocate_2c(zSq,M0,M0,infoloc)
    call wallocate_2c(work1,N,M0,infoloc)
    call wallocate_2c(work2,N,M0,infoloc)
    call wallocate_2c(Az,band,N,infoloc)
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
                s=min(klz,N-i)
                call CCOPY(s,A(2,i),1,Az(klz+2,i),1)
                call CSCAL(s,-ONEC,Az(klz+2,i),1)
                call CCOPY(s,Az(klz+2,i),1,Az(klz,i+1),2*klz)
                call CLACGV(s, Az(klz,i+1), 2*klz )
                Az(klz+1,i)=-A(1,i)+Ze
             enddo
             Az(klz+1,N)=-A(1,N)+Ze

          elseif ((UPLO=='U').or.(UPLO=='u')) then

             do i=1,N-1
                s=min(klz,N-i)
                call CCOPY(s,A(klz,i+1),LDA-1,Az(klz,i+1),2*klz)
                call CSCAL(s,-ONEC,Az(klz,i+1),2*klz)
                call CCOPY(s,Az(klz,i+1),2*klz,Az(klz+2,i),1)
                call CLACGV(s, Az(klz+2,i),1)
                Az(klz+1,i)=-A(klz+1,i)+Ze
             enddo
             Az(klz+1,N)=-A(klz+1,N)+Ze


          else !UPLO=F

             do i=1,N
                call CCOPY(band,A(1,i),1,Az(1,i),1)
                call CSCAL(band,-ONEC,Az(1,i),1)
                Az(klz+1,i)=Az(klz+1,i)+Ze
             end do


          end if



          nzero=SZERO
          call CGBALU(N,klz,klz,Az,band,nzero,norm,infoloc)
          if (infoloc/=0) then
             info=-2
             return
          end if


       case(11) !!solve the linear system (ZeB-A)x=work2(1:N,1:M0) result in to work2



          call CTBSM('L','N','U',N,M0,klz,Az(klz+1,1),band,work2,N)
          call CTBSM('U','N','N',N,M0,klz,Az(1,1),band,work2,N)




       case(21) !!solve the linear system (ZeB-A)^H x=work2(1:N,1:M0) result in to work2

          call CTBSM('U','C','N',N,M0,klz,Az(1,1),band,work2,N)
          call CTBSM('L','C','U',N,M0,klz,Az(klz+1,1),band,work2,N)


       case(30) !! perform multiplication A*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call CHBMM(UPLO,n,fpm(25),kla,ONEC,A(1,1),LDA,X(1,fpm(24)),N,ZEROC,work1(1,fpm(24)),N)

       case(40) !! perform multiplication B*x(1:N,fpm(24):fpm(24)+fpm(25)-1) result in work(1:N,fpm(24)+fpm(25)-1)

          call CLACPY( 'F', N, fpm(25),X(1,fpm(24)) , N, work1(1,fpm(24)), N )


       end select
    end do



    call wdeallocate_2c(zAq)
    call wdeallocate_2c(zSq)
    call wdeallocate_2c(work2)
    call wdeallocate_2c(Az)
    call wdeallocate_2c(work1)

  end subroutine cfeast_hbev


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
















