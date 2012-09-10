module basis_aux
use types, only: dp
use qc, only: coulomb_repulsion
use scf, only: ijkl2intindex
implicit none
private
public get_2ints_size, ints_shell_set

contains

subroutine ints_shell_set(A, B, C, D, &
        powera, powerb, powerc, powerd, &
        nprima, nprimb, nprimc, nprimd, &
        coefa, coefb, coefc, coefd, &
        alphaa, alphab, alphac, alphad, &
        ishell, jshell, kshell, lshell, &
        ilambda, jlambda, klambda, llambda, &
        r)
real(dp), intent(in), dimension(3) :: A, B, C, D
integer, intent(in), dimension(:, :) :: powera, powerb, powerc, powerd
integer, intent(in) :: nprima, nprimb, nprimc, nprimd
real(dp), intent(in) :: alphaa(:), alphab(:), alphac(:), alphad(:)
real(dp), intent(in) :: coefa(:), coefb(:), coefc(:), coefd(:)
integer, intent(in) :: ishell, jshell, kshell, lshell
integer, intent(in) :: ilambda, jlambda, klambda, llambda
real(dp), intent(inout) :: r(:)

integer :: ip, jp, kp, lp
integer :: i, j, k, l
real(dp) :: integ, incr
logical :: can_compute
!can_compute = (lshell(is) >= lshell(js) .and. &
!        lshell(ks) >= lshell(ls) &
!    .and. lshell(is)+lshell(js) <= lshell(ks)+lshell(ls))
!if (.not. can_compute) cycle

do i = ishell, ishell + (ilambda+1)*(ilambda+2)/2-1
do j = jshell, jshell + (jlambda+1)*(jlambda+2)/2-1
do k = kshell, kshell + (klambda+1)*(klambda+2)/2-1
do l = lshell, lshell + (llambda+1)*(llambda+2)/2-1

                    can_compute = (i >= j .and. k >= l .and. &
                        i*(i-1)/2+j >= k*(k-1)/2+l)
                    if (.not. can_compute) cycle
integ = 0
! Loop over primitive functions
do ip = 1, nprima
    do jp = 1, nprimb
        do kp = 1, nprimc
            do lp = 1, nprimd
                incr = coulomb_repulsion( &
                A(1), A(2), A(3), powera(1, i), powera(2, i), powera(3, i), alphaa(ip), &
                B(1), B(2), B(3), powerb(1, j), powerb(2, j), powerb(3, j), alphab(jp), &
                C(1), C(2), C(3), powerc(1, k), powerc(2, k), powerc(3, k), alphac(kp), &
                D(1), D(2), D(3), powerd(1, l), powerd(2, l), powerd(3, l), alphad(lp))
                integ = integ + &
                    coefa(ip+nprima*(i-ishell)) * &
                    coefb(jp+nprimb*(j-jshell)) * &
                    coefc(kp+nprimc*(k-kshell)) * &
                    coefd(lp+nprimd*(l-lshell)) * &
                    incr
            end do
        end do
    end do
end do

r(ijkl2intindex(i,j,k,l)) = integ
end do
end do
end do
end do
end subroutine

integer pure function get_2ints_size(n) result(r)
integer, intent(in) :: n
integer :: m
! Number of ij combinations
m = n*(n+1)/2
! Number of ij vs kl combinations in (ij|kl)
r = m*(m+1)/2
end function

end module
