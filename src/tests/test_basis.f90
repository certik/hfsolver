program test_basis
use basis, only: get_cartesian_shell, getstart, get_2ints_size
use utils, only: assert
implicit none

integer :: l, m, n
integer, allocatable, dimension(:) :: lpower, mpower, npower
l = 0
allocate(lpower((l+1)*(l+2)/2), mpower((l+1)*(l+2)/2), npower((l+1)*(l+2)/2))
call get_cartesian_shell(l, lpower, mpower, npower)
call assert(all(lpower == [0]))
call assert(all(mpower == [0]))
call assert(all(npower == [0]))
deallocate(lpower, mpower, npower)

l = 1
allocate(lpower((l+1)*(l+2)/2), mpower((l+1)*(l+2)/2), npower((l+1)*(l+2)/2))
call get_cartesian_shell(l, lpower, mpower, npower)
call assert(all(lpower == [1, 0, 0]))
call assert(all(mpower == [0, 1, 0]))
call assert(all(npower == [0, 0, 1]))
deallocate(lpower, mpower, npower)

l = 2
allocate(lpower((l+1)*(l+2)/2), mpower((l+1)*(l+2)/2), npower((l+1)*(l+2)/2))
call get_cartesian_shell(l, lpower, mpower, npower)
call assert(all(lpower == [2, 1, 1, 0, 0, 0]))
call assert(all(mpower == [0, 1, 0, 2, 1, 0]))
call assert(all(npower == [0, 0, 1, 0, 1, 2]))
deallocate(lpower, mpower, npower)

l = 3
allocate(lpower((l+1)*(l+2)/2), mpower((l+1)*(l+2)/2), npower((l+1)*(l+2)/2))
call get_cartesian_shell(l, lpower, mpower, npower)
call assert(all(lpower == [3, 2, 2, 1, 1, 1, 0, 0, 0, 0]))
call assert(all(mpower == [0, 1, 0, 2, 1, 0, 3, 2, 1, 0]))
call assert(all(npower == [0, 0, 1, 0, 1, 2, 0, 1, 2, 3]))
deallocate(lpower, mpower, npower)

l = 4
allocate(lpower((l+1)*(l+2)/2), mpower((l+1)*(l+2)/2), npower((l+1)*(l+2)/2))
call get_cartesian_shell(l, lpower, mpower, npower)
call assert(all(lpower == [4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0]))
call assert(all(mpower == [0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0]))
call assert(all(npower == [0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4]))
deallocate(lpower, mpower, npower)

! -------------------------------------------------------
call assert(all(getstart([6, 3, 3, 1, 1, 1]) == [1, 7, 10, 13, 14, 15]))
call assert(all(getstart([6, 6, 6, 3, 3, 1, 1, 1]) == &
    [1, 7, 13, 19, 22, 25, 26, 27]))

! -------------------------------------------------------
n = 15
m = n*16/2
call assert(get_2ints_size(n) == m*(m+1)/2)

end program
