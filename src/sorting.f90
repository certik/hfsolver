module sorting

! Module for sorting arrays.
! Based on code written by John E. Pask, LLNL.

use types, only: dp
use utils, only: stop_error
implicit none
private
public sort, sortpairs, argsort, iargsort_quicksort

! overload argsort
interface argsort
    module procedure iargsort, rargsort
end interface

! overload sort
interface sort
    module procedure sortNums, sortINums, sortVecs
end interface

! overload sortpairs
interface sortpairs
    module procedure sortNumNumPairs, sortINumCNumPairs, &
                     sortNumVecPairs, sortINumVecPairs, sortNumCVecPairs, &
                     sortNumMatPairs, sortINumMatPairs
end interface

contains

subroutine sortNums(nums)
! sorts array of numbers, nums, from smallest to largest
real(dp), intent(inout):: nums(:)   ! array of numbers
nums = nums(argsort(nums))
end subroutine

subroutine sortINums(nums)
! sorts array of inegers, nums, from smallest to largest
integer, intent(inout):: nums(:)    ! array of numbers
nums = nums(argsort(nums))
end subroutine

subroutine sortVecs(vecs)
! sorts array of vectors, vecs, by length, from smallest to largest
real(dp), intent(inout):: vecs(:,:) ! array of vectors: vecs(i,j) = ith comp. of jth vec.
real(dp) len2(size(vecs,2))         ! array of squares of vector lengths
integer i
do i=1,size(len2)
   len2(i)=dot_product(vecs(:,i),vecs(:,i))
end do
call sortpairs(len2,vecs)
end subroutine

subroutine sortNumNumPairs(nums1, nums2)
! sorts arrays of numbers, nums1 and nums2, according to increasing nums1
real(dp), intent(inout):: nums1(:), nums2(:)  ! arrays of numbers
integer :: a(size(nums1))
if (size(nums1) /= size(nums2)) then
    call stop_error("SORTPAIRS ERROR: arrays must be of same length.")
end if
a = argsort(nums1)
nums1 = nums1(a)
nums2 = nums2(a)
end subroutine

subroutine sortINumCNumPairs(nums1, nums2)
! sorts arrays of numbers, nums1 and nums2, according to increasing nums1
integer, intent(inout):: nums1(:)            ! array of integers
complex(dp), intent(inout):: nums2(:)        ! array of complex numbers
integer :: a(size(nums1))
if (size(nums1) /= size(nums2)) then
    call stop_error("SORTPAIRS ERROR: arrays must be of same length.")
end if
a = argsort(nums1)
nums1 = nums1(a)
nums2 = nums2(a)
end subroutine

subroutine sortNumVecPairs(nums, vecs)
! sorts arrays of numbers, nums, and vectors, vecs, according to increasing nums
real(dp), intent(inout):: nums(:)   ! array of numbers
real(dp), intent(inout):: vecs(:,:) ! array of vectors: vecs(i,j) = ith comp. of jth vec.
integer :: a(size(nums))
if (size(nums) /= size(vecs,2)) then
    call stop_error("SORTPAIRS ERROR: arrays must be of same length.")
end if
a = argsort(nums)
nums = nums(a)
vecs = vecs(:, a)
end subroutine

subroutine sortINumVecPairs(nums, vecs)
! sorts arrays of integers, nums, and vectors, vecs, according to increasing nums
integer, intent(inout):: nums(:)    ! array of numbers
real(dp), intent(inout):: vecs(:,:) ! array of vectors: vecs(i,j) = ith comp. of jth vec.
integer :: a(size(nums))
if (size(nums) /= size(vecs,2)) then
    call stop_error("SORTPAIRS ERROR: arrays must be of same length.")
end if
a = argsort(nums)
nums = nums(a)
vecs = vecs(:, a)
end subroutine

subroutine sortNumCVecPairs(nums, vecs)
! sorts arrays of numbers, nums, and complex vectors, vecs, according to increasing nums
real(dp), intent(inout):: nums(:)   ! array of numbers
complex(dp), intent(inout):: vecs(:,:) ! array of vectors: vecs(i,j) = ith comp. of jth vec.
integer :: a(size(nums))
if (size(nums) /= size(vecs,2)) then
    call stop_error("SORTPAIRS ERROR: arrays must be of same length.")
end if
a = argsort(nums)
nums = nums(a)
vecs = vecs(:, a)
end subroutine

subroutine sortNumMatPairs(nums, mats)
! sorts arrays of numbers, nums, and matrices, mats, according to increasing nums
real(dp), intent(inout):: nums(:)         ! array of numbers
real(dp), intent(inout):: mats(:,:,:)     ! array of matrices: mats(i,j,n) = (i,j) comp. of nth mat.
integer :: a(size(nums))
if (size(nums) /= size(mats,3)) then
    call stop_error("SORTPAIRS ERROR: arrays must be of same length.")
end if
a = argsort(nums)
nums = nums(a)
mats = mats(:, :, a)
end subroutine

subroutine sortINumMatPairs(nums,mats)
! sorts arrays of integers, nums, and matrices, mats, according to increasing nums
integer, intent(inout):: nums(:)          ! array of numbers
real(dp), intent(inout):: mats(:,:,:)     ! array of matrices: mats(i,j,n) = (i,j) comp. of nth mat.
integer :: a(size(nums))
if (size(nums) /= size(mats,3)) then
    call stop_error("SORTPAIRS ERROR: arrays must be of same length.")
end if
a = argsort(nums)
nums = nums(a)
mats = mats(:, :, a)
end subroutine

pure elemental subroutine swap_int(x,y)
  integer, intent(in out) :: x,y
  integer :: z
  z = x
  x = y
  y = z
end subroutine

pure subroutine interchange_sort_map_int(vec,map)
  integer, intent(in out) :: vec(:)
  integer, intent(in out) :: map(:)
  integer :: i,j
  do i = 1,size(vec) - 1
     j = minloc(vec(i:),1)
     if (j > 1) then
        call swap_int(vec(i),vec(i + j - 1))
        call swap_int(map(i),map(i + j - 1))
     end if
  end do
end subroutine

pure function iargsort_quicksort(vec_) result(map)
  integer, intent(in) :: vec_(:)
  integer :: map(size(vec_))
  integer, parameter :: levels = 300
  integer, parameter :: max_interchange_sort_size = 20
  integer :: i,left,right,l_bound(levels),u_bound(levels)
  integer :: pivot
  integer :: vec(size(vec_))

  vec = vec_

  forall(i=1:size(vec)) map(i) = i

  l_bound(1) = 1
  u_bound(1) = size(vec)
  i = 1
  do while(i >= 1)
     left = l_bound(i)
     right = u_bound(i)
     if (right - left < max_interchange_sort_size) then
        if (left < right) call interchange_sort_map_int(vec(left:right),map(left:right))
        i = i - 1
     else
        pivot = (vec(left) + vec(right)) / 2
        left = left - 1
        right = right + 1
        do
           do
              left = left + 1
              if (vec(left) >= pivot) exit
           end do
           do
              right = right - 1
              if (vec(right) <= pivot) exit
           end do
           if (left < right) then
              call swap_int(vec(left),vec(right))
              call swap_int(map(left),map(right))
           elseif(left == right) then
              if (left == l_bound(i)) then
                 left = left + 1
              else
                 right = right - 1
              end if
              exit
           else
              exit
           end if
           end do
           u_bound(i + 1) = u_bound(i)
           l_bound(i + 1) = left
           u_bound(i) = right
           i = i + 1
     end if
  end do
end function

function iargsort(a) result(b)
! Returns the indices that would sort an array.
!
! Arguments
! ---------
!
integer, intent(in):: a(:)    ! array of numbers
integer :: b(size(a))         ! indices into the array 'a' that sort it
!
! Example
! -------
!
! iargsort([10, 9, 8, 7, 6])   ! Returns [5, 4, 3, 2, 1]

integer :: N                           ! number of numbers/vectors
integer :: i,imin,relimin(1)           ! indices: i, i of smallest, relative imin
integer :: temp                        ! temporary
integer :: a2(size(a))
a2 = a
N=size(a)
do i = 1, N
    b(i) = i
end do
do i = 1, N-1
    ! find ith smallest in 'a'
    relimin = minloc(a2(i:))
    imin = relimin(1) + i - 1
    ! swap to position i in 'a' and 'b', if not already there
    if (imin /= i) then
        temp = a2(i); a2(i) = a2(imin); a2(imin) = temp
        temp = b(i); b(i) = b(imin); b(imin) = temp
    end if
end do
end function

function rargsort(a) result(b)
! Returns the indices that would sort an array.
!
! Arguments
! ---------
!
real(dp), intent(in):: a(:)   ! array of numbers
integer :: b(size(a))         ! indices into the array 'a' that sort it
!
! Example
! -------
!
! rargsort([4.1_dp, 2.1_dp, 2.05_dp, -1.5_dp, 4.2_dp]) ! Returns [4, 3, 2, 1, 5]

integer :: N                           ! number of numbers/vectors
integer :: i,imin,relimin(1)           ! indices: i, i of smallest, relative imin
integer :: temp1                        ! temporary
real(dp) :: temp2
real(dp) :: a2(size(a))
a2 = a
N=size(a)
do i = 1, N
    b(i) = i
end do
do i = 1, N-1
    ! find ith smallest in 'a'
    relimin = minloc(a2(i:))
    imin = relimin(1) + i - 1
    ! swap to position i in 'a' and 'b', if not already there
    if (imin /= i) then
        temp2 = a2(i); a2(i) = a2(imin); a2(imin) = temp2
        temp1 = b(i); b(i) = b(imin); b(imin) = temp1
    end if
end do
end function

end module
