module erd
implicit none

interface
    subroutine  ERD__GENER_ERI_BATCH (IMAX,ZMAX, &
            NALPHA,NCOEFF,NCSUM, &
            NCGTO1,NCGTO2,NCGTO3,NCGTO4, &
            NPGTO1,NPGTO2,NPGTO3,NPGTO4, &
            SHELL1,SHELL2,SHELL3,SHELL4, &
            X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4, &
            ALPHA,CC,CCBEG,CCEND, &
            SPHERIC, &
            SCREEN, &
            ICORE, &
                NBATCH, &
                NFIRST, &
                ZCORE )
    implicit none
    LOGICAL     SCREEN
    LOGICAL     SPHERIC

    INTEGER     IMAX,ZMAX
    INTEGER     NALPHA,NCOEFF,NCSUM
    INTEGER     NBATCH,NFIRST
    INTEGER     NCGTO1,NCGTO2,NCGTO3,NCGTO4
    INTEGER     NPGTO1,NPGTO2,NPGTO3,NPGTO4
    INTEGER     SHELL1,SHELL2,SHELL3,SHELL4

    INTEGER     CCBEG (1:NCSUM)
    INTEGER     CCEND (1:NCSUM)
    INTEGER     ICORE (1:IMAX)

    DOUBLE PRECISION  X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4

    DOUBLE PRECISION  ALPHA (1:NALPHA)
    DOUBLE PRECISION  CC    (1:NCOEFF)
    DOUBLE PRECISION  ZCORE (1:ZMAX)
    end subroutine
end interface

end module

program test_erd
use erd, only: erd__gener_eri_batch
implicit none

integer, parameter :: dp = kind(0.d0)

integer, parameter :: IMAX=100, ZMAX=10000
integer, parameter :: NALPHA=4, NCOEFF=4, NCSUM=4

LOGICAL     SCREEN
LOGICAL     SPHERIC
INTEGER     NCGTO1,NCGTO2,NCGTO3,NCGTO4
INTEGER     NPGTO1,NPGTO2,NPGTO3,NPGTO4
INTEGER     SHELL1,SHELL2,SHELL3,SHELL4
INTEGER     CCBEG (1:NCSUM)
INTEGER     CCEND (1:NCSUM)
INTEGER     ICORE (1:IMAX)
real(dp) X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4
real(dp) ALPHA(1:NALPHA)
real(dp) CC(1:NCOEFF)

integer NBATCH,NFIRST
real(dp) ZCORE(1:ZMAX)

print *, "Start"
screen = .false.
spheric = .false.
x1 = 0; y1 = 0; z1 = 0
x2 = 0; y2 = 0; z2 = 0
x3 = 0; y3 = 0; z3 = 0
x4 = 0; y4 = 0; z4 = 0
ncgto1 = 1; ncgto2 = 1; ncgto3 = 1; ncgto4 = 1
npgto1 = 1; npgto2 = 1; npgto3 = 1; npgto4 = 1
shell1 = 1; shell2 = 1; shell3 = 1; shell4 = 1
alpha = [1._dp, 2._dp, 2._dp, 2._dp]
cc = [1, 1, 1, 1]
ccbeg = [1, 1, 1, 1]
ccend = [1, 1, 1, 1]
call ERD__GENER_ERI_BATCH (IMAX,ZMAX, &
            NALPHA,NCOEFF,NCSUM, &
            NCGTO1,NCGTO2,NCGTO3,NCGTO4, &
            NPGTO1,NPGTO2,NPGTO3,NPGTO4, &
            SHELL1,SHELL2,SHELL3,SHELL4, &
            X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4, &
            ALPHA,CC,CCBEG,CCEND, &
            SPHERIC, &
            SCREEN, &
            ICORE, &
                NBATCH, &
                NFIRST, &
                ZCORE )
print *, "Done"
print *, "nbatch =", NBATCH
print *, "nfirst =", NFIRST
! Integrals:
!print *, zcore(nfirst:nfirst+nbatch)
end
