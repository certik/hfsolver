program test_ijkl
use scf, only: ijkl2intindex, ijkl2intindex2
use utils, only: assert
implicit none

call testijkl(1)
call testijkl(2)
call testijkl(10)
call testijkl(50)

call testijkl2(1)
call testijkl2(2)
call testijkl2(9)
call testijkl2(15)
! Too slow:
!call testijkl2(50)

contains

subroutine testijkl(n)
integer, intent(in) :: n
integer :: i, j, k, l, ijkl
ijkl = 1
do i = 1, n
    do j = 1, i
        do k = 1, n
            do l = 1, k
                if ((i-1)*i/2 + j < (k-1)*k/2 + l) cycle
                call assert(ijkl2intindex(i, j, k, l) == ijkl)
                call assert(ijkl2intindex(j, i, k, l) == ijkl)
                call assert(ijkl2intindex(i, j, l, k) == ijkl)
                call assert(ijkl2intindex(j, i, l, k) == ijkl)
                call assert(ijkl2intindex(k, l, i, j) == ijkl)
                call assert(ijkl2intindex(l, k, i, j) == ijkl)
                call assert(ijkl2intindex(k, l, j, i) == ijkl)
                call assert(ijkl2intindex(l, k, j, i) == ijkl)
                ijkl = ijkl + 1
            end do
        end do
    end do
end do
ijkl = ijkl - 1
! Test the total number of matrix elements:
call assert(ijkl == n * (n**3 + 2*n**2 + 3*n + 2) / 8)
end subroutine


subroutine testijkl2(n)
integer, intent(in) :: n
integer :: i, j, k, l, ijkl
ijkl = 1
do i = 1, n
    do j = 1, i
        do k = 1, i
            do l = 1, i
                if (i == j .and. k < l) cycle
                if (i == k .and. j < l) cycle
                if (i == l .and. j < k) cycle
                call assert(ijkl2intindex2(i, j, k, l) == ijkl)
                call assert(ijkl2intindex2(j, i, l, k) == ijkl)
                call assert(ijkl2intindex2(k, l, i, j) == ijkl)
                call assert(ijkl2intindex2(l, k, j, i) == ijkl)
                ijkl = ijkl + 1
            end do
        end do
    end do
end do
ijkl = ijkl - 1
! Test the total number of matrix elements:
call assert(ijkl == n**2 * (n**2 + 3) / 4)
end subroutine

end program
