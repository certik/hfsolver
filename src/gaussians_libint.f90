module gaussians

! Two particle Gaussian integrals are calculated using libint

use types, only: dp
use basis_aux, only: get_2ints_size, ints_shell_set
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
                if (ilambda == 0 .and. jlambda == 0 .and. klambda == 0 .and. &
                    llambda == 0) then
                    call cpu_time(t1)
                    call ints_shell_set( &
                        center(:, i), center(:, j), center(:, k), center(:, l),&
                        power, power, power, power, &
                        nprim(i), nprim(j), nprim(k), nprim(l), &
                        coef(istart(i):istart(i)+nprim(i)*(ilambda+1)*(ilambda+2)/2-1), &
                        coef(istart(j):istart(j)+nprim(j)*(jlambda+1)*(jlambda+2)/2-1), &
                        coef(istart(k):istart(k)+nprim(k)*(klambda+1)*(klambda+2)/2-1), &
                        coef(istart(l):istart(l)+nprim(l)*(llambda+1)*(llambda+2)/2-1), &
                        alpha(istart(i):istart(i)+nprim(i)-1), &
                        alpha(istart(j):istart(j)+nprim(j)-1), &
                        alpha(istart(k):istart(k)+nprim(k)-1), &
                        alpha(istart(l):istart(l)+nprim(l)-1), &
                        i, j, k, l, &
                        ilambda, jlambda, klambda, llambda, &
                        int2)
                    call cpu_time(t2)
                    tS = tS + t2-t1
                else
                    call cpu_time(t1)
                    call libint_getInts2( &
                        center(:, i), center(:, j), center(:, k), center(:, l),&
                        nprim(i), nprim(j), nprim(k), nprim(l), &
                        coef(istart(i):istart(i)+nprim(i)*(ilambda+1)*(ilambda+2)/2-1), &
                        coef(istart(j):istart(j)+nprim(j)*(jlambda+1)*(jlambda+2)/2-1), &
                        coef(istart(k):istart(k)+nprim(k)*(klambda+1)*(klambda+2)/2-1), &
                        coef(istart(l):istart(l)+nprim(l)*(llambda+1)*(llambda+2)/2-1), &
                        alpha(istart(i):istart(i)+nprim(i)-1), &
                        alpha(istart(j):istart(j)+nprim(j)-1), &
                        alpha(istart(k):istart(k)+nprim(k)-1), &
                        alpha(istart(l):istart(l)+nprim(l)-1), &
                        i, j, k, l, &
                        ilambda, jlambda, klambda, llambda, &
                        0, &
                        size(int2), int2)
                    call cpu_time(t2)
                    tP = tP + t2-t1
                end if
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


end module
