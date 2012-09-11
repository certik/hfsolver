module gaussians

! Two particle Gaussian integrals are calculated using libint

use types, only: dp
use utils, only: stop_error
use scf, only: ijkl2intindex
implicit none
private
public getints2

interface

    subroutine libint_getInts2(A, B, C, D, &
            nprima, nprimb, nprimc, nprimd, &
            coefa, coefb, coefc, coefd, &
            alphaa, alphab, alphac, alphad, &
            ishell, jshell, kshell, lshell, &
            ilambda, jlambda, klambda, llambda, &
            permut, &
            n2, r) bind(c, name="getInts2")
    use iso_c_binding, only: c_int, c_double
    implicit none
    real(c_double), intent(in), dimension(3) :: A, B, C, D
    integer(c_int), value, intent(in) :: nprima, nprimb, nprimc, nprimd
    real(c_double), intent(in) :: alphaa(nprima), alphab(nprimb), &
        alphac(nprimc), alphad(nprimd)
    integer(c_int), value, intent(in) :: ilambda, jlambda, klambda, llambda
    real(c_double), intent(in) :: &
        coefa(nprima*(ilambda+1)*(ilambda+2)/2), &
        coefb(nprimb*(jlambda+1)*(jlambda+2)/2), &
        coefc(nprimc*(klambda+1)*(klambda+2)/2), &
        coefd(nprimd*(llambda+1)*(llambda+2)/2)
    integer(c_int), value, intent(in) :: ishell, jshell, kshell, lshell
    integer(c_int), value, intent(in) :: permut
    integer(c_int), value, intent(in) :: n2
    real(c_double), intent(inout) :: r(n2)
    end subroutine

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

contains

subroutine getints2(nprim, istart, center, power, coef, alpha, int2)
integer, intent(in) :: nprim(:), istart(:)
real(dp), intent(in) :: center(:, :)
integer, intent(in) :: power(:, :)
real(dp), intent(in) :: coef(:), alpha(:)
real(dp), intent(out) :: int2(:)

integer :: i, j, k, l ! these go over shells, not basis functions
integer :: ilambda
integer :: jlambda
integer :: klambda
integer :: llambda
real(dp) :: tS, tP, t1, t2, tb, te
int2 = 0
call cpu_time(tb)
tS = 0
tP = 0
i = 1
do
    ilambda = sum(power(:, i))
    j = 1
    do
        jlambda = sum(power(:, j))
        k = 1
        do
            klambda = sum(power(:, k))
            l = 1
            do
                llambda = sum(power(:, l))
                call cpu_time(t1)
                call do_erd( &
                    center(:, i), center(:, j), center(:, k), center(:, l),&
                    nprim(i), nprim(j), nprim(k), nprim(l), &
                    coef(istart(i):istart(i)+nprim(i)-1), &
                    coef(istart(j):istart(j)+nprim(j)-1), &
                    coef(istart(k):istart(k)+nprim(k)-1), &
                    coef(istart(l):istart(l)+nprim(l)-1), &
                    alpha(istart(i):istart(i)+nprim(i)-1), &
                    alpha(istart(j):istart(j)+nprim(j)-1), &
                    alpha(istart(k):istart(k)+nprim(k)-1), &
                    alpha(istart(l):istart(l)+nprim(l)-1), &
                    i, j, k, l, &
                    ilambda, jlambda, klambda, llambda, &
                    int2)
                call cpu_time(t2)
                tP = tP + t2-t1
                l = l + (llambda+1)*(llambda+2)/2
                if (l > size(nprim)) exit
            end do ! l
            k = k + (klambda+1)*(klambda+2)/2
            if (k > size(nprim)) exit
        end do ! k
        j = j + (jlambda+1)*(jlambda+2)/2
        if (j > size(nprim)) exit
    end do ! j
    i = i + (ilambda+1)*(ilambda+2)/2
    if (i > size(nprim)) exit
end do ! i
call cpu_time(te)
print *, "TIMING:"
print "(a,f10.6)", "    (ss|ss)    : ", tS
print "(a,f10.6)", "    (..|..)    : ", tP
print "(a,f10.6)", "    sum        : ", tS+tP
print "(a,f10.6)", "  direct total : ", te-tb
end subroutine

subroutine do_erd(A, B, C, D, &
        nprima, nprimb, nprimc, nprimd, &
        coefa, coefb, coefc, coefd, &
        alphaa, alphab, alphac, alphad, &
        ishell, jshell, kshell, lshell, &
        ilambda, jlambda, klambda, llambda, &
        r)
real(dp), intent(in), dimension(3) :: A, B, C, D
integer, intent(in) :: nprima, nprimb, nprimc, nprimd
real(dp), intent(in) :: alphaa(:), alphab(:), alphac(:), alphad(:)
real(dp), intent(in) :: coefa(:), coefb(:), coefc(:), coefd(:)
integer, intent(in) :: ishell, jshell, kshell, lshell
integer, intent(in) :: ilambda, jlambda, klambda, llambda
real(dp), intent(inout) :: r(:)

integer, parameter :: IMAX=100000, ZMAX=100000
integer :: icore(IMAX)
real(dp) :: zcore(ZMAX)

logical :: screen, spheric
integer :: ncgto1,ncgto2,ncgto3,ncgto4
integer :: npgto1,npgto2,npgto3,npgto4
integer :: shell1,shell2,shell3,shell4
! We only use one contraction per shell so far:
integer, parameter :: NCSUM=4 ! Number of contractions (sum of all shells)
integer :: ccbeg(NCSUM), ccend(NCSUM)
real(dp) :: X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4
real(dp), dimension(nprima + nprimb + nprimc + nprimd) :: alpha, cc
integer :: nbatch, nfirst
integer :: i, j, k, l, counter
logical :: can_compute

print *, "Start"
print *, "shell:", ishell, jshell, kshell, lshell
print *, "lambda:", ilambda, jlambda, klambda, llambda
print *, "prim:", nprima, nprimb, nprimc, nprimd
print *, coefa
print *, coefb
print *, coefc
print *, coefd
screen = .false.
spheric = .false.
x1 = A(1); y1 = A(2); z1 = A(3)
x2 = B(1); y2 = B(2); z2 = B(3)
x3 = C(1); y3 = C(2); z3 = C(3)
x4 = D(1); y4 = D(2); z4 = D(3)
! We only use one contraction per shell so far:
ncgto1 = 1; ncgto2 = 1; ncgto3 = 1; ncgto4 = 1
npgto1 = nprima; npgto2 = nprimb; npgto3 = nprimc; npgto4 = nprimd
shell1 = ilambda; shell2 = jlambda; shell3 = klambda; shell4 = llambda
! Pack the alpha and coef coefficients:
alpha = [alphaa, alphab, alphac, alphad]
cc = [coefa, coefb, coefc, coefd]
! For non-segmented contractions, this should be correct:
ccbeg = [1, 1, 1, 1]
ccend = [npgto1, npgto2, npgto3, npgto4]
call ERD__GENER_ERI_BATCH(size(icore), size(zcore), &
            size(alpha),size(cc),size(ccbeg), &
            NCGTO1,NCGTO2,NCGTO3,NCGTO4, &
            NPGTO1,NPGTO2,NPGTO3,NPGTO4, &
            SHELL1,SHELL2,SHELL3,SHELL4, &
            X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4, &
            alpha, cc, ccbeg, ccend, &
            spheric, screen, icore, &
                nbatch, nfirst, zcore)
print *, "Done"
print *, "nbatch =", NBATCH
print *, "nfirst =", NFIRST
print *, "size(r) =", size(r)
if (nbatch == 0) then
    ! All integrals are zero
    return
end if
if (nbatch /=  (ilambda+1)*(ilambda+2)/2 * &
        (jlambda+1)*(jlambda+2)/2 * &
        (klambda+1)*(klambda+2)/2 * &
        (llambda+1)*(llambda+2)/2) then
    call stop_error("Incorrect number of integrals returned")
end if
counter = 0
do i = ishell, ishell + (ilambda+1)*(ilambda+2)/2-1
do j = jshell, jshell + (jlambda+1)*(jlambda+2)/2-1
do k = kshell, kshell + (klambda+1)*(klambda+2)/2-1
do l = lshell, lshell + (llambda+1)*(llambda+2)/2-1
    can_compute = (i >= j .and. k >= l .and. &
        i*(i-1)/2+j >= k*(k-1)/2+l)
! TODO: handle this somehow:
!    if (.not. can_compute) cycle

    print *, "setting:", i, j, k, l, ijkl2intindex(i, j, k, l), zcore(nfirst+counter)
    r(ijkl2intindex(i,j,k,l)) = zcore(nfirst+counter)
    counter = counter + 1
end do
end do
end do
end do
end subroutine

end module
