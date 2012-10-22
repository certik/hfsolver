module feast
implicit none

integer, parameter:: dp=kind(0.d0)

interface

    subroutine feastinit(fpm)
    integer,dimension(*) :: fpm
    end subroutine

    subroutine dfeast_syev(UPLO,N,A,LDA,fpm,epsout,loop,Emin,Emax,M0, &
            E,X,mode,res,info)
    implicit none
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
    end subroutine

    subroutine dfeast_sygv(UPLO,N,A,LDA,B,LDB,fpm,epsout,loop,Emin,Emax,M0, &
            E,X,mode,res,info)
    implicit none
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
    end subroutine

    subroutine dfeast_srci(ijob,N,Ze,work,workc,Aq,Sq,fpm,epsout,loop, &
            Emin,Emax,M0,lambda,q,mode,res,info)
    implicit none
    integer :: ijob,N,M0
    complex(kind=(kind(1.0d0))) :: Ze
    double precision, dimension(N,*) ::work
    complex(kind=(kind(1.0d0))),dimension(N,*):: workc
    double precision,dimension(M0,*):: Aq,Sq
    integer,dimension(*) :: fpm
    double precision :: epsout
    integer :: loop
    double precision :: Emin,Emax
    double precision,dimension(*)  :: lambda
    double precision,dimension(N,*):: q
    integer :: mode
    double precision,dimension(*) :: res
    integer :: info
    end subroutine

end interface

contains

end module
