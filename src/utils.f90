module utils

! Various general utilities.
! Based on a code by John E. Pask, LLNL.

use types, only: dp
implicit none
private
public upcase, lowcase, whitechar, blank, num_strings, getstring, &
    stop_error, arange, loadtxt, savetxt, newunit, assert, str, init_random, &
    zeros, mesh_exp, linspace, clock, strfmt, get_int_arg, get_float_arg, &
    allocate_mold

interface str
    module procedure str_int, str_real, str_real_n
end interface

interface strfmt
    module procedure strfmt_int, strfmt_real
end interface

interface allocate_mold
    module procedure allocate_mold_real3d
    module procedure allocate_mold_complex3d
end interface

contains

function upcase(s) result(t)
! Returns string 's' in uppercase
character(*), intent(in) :: s
character(len(s)) :: t
integer :: i, diff
t = s; diff = ichar('A')-ichar('a')
do i = 1, len(t)
    if (ichar(t(i:i)) >= ichar('a') .and. ichar(t(i:i)) <= ichar('z')) then
        ! if lowercase, make uppercase
        t(i:i) = char(ichar(t(i:i)) + diff)
    end if
end do
end function

function lowcase(s) result(t)
! Returns string 's' in lowercase
character(*), intent(in) :: s
character(len(s)) :: t
integer :: i, diff
t = s; diff = ichar('A')-ichar('a')
do i = 1, len(t)
    if (ichar(t(i:i)) >= ichar('A') .and. ichar(t(i:i)) <= ichar('Z')) then
        ! if uppercase, make lowercase
        t(i:i) = char(ichar(t(i:i)) - diff)
    end if
end do
end function

logical function whitechar(char) ! white character
! returns .true. if char is space (32) or tab (9), .false. otherwise
character, intent(in) :: char
if (iachar(char) == 32 .or. iachar(char) == 9) then
    whitechar = .true.
else
    whitechar = .false.
end if
end function

logical function blank(string)
! Returns true if string contains only white characters
character(*), intent(in) :: string
integer :: i
do i = 1, len(string)
    if (.not. whitechar(string(i:i))) exit
end do
blank = (i>len(string))
end function

integer function num_strings(s) result(n)
! Returns number of substrings contained in input string 's' delimited
! by white space.
character(*), intent(in) :: s    ! input string
character(len(s)+2) :: t         ! temporary string to facilitate analysis
integer :: i
t = " " // s // " "
n = 0
do i = 1, len(t)-1
    if (whitechar(t(i:i)) .and. .not. whitechar(t(i+1:i+1))) n = n + 1
end do
end function

!--------------------------------------------------------------------------------------------------!

subroutine getstring(s,is,ss)
! Returns first substring ss in string s, delimited by white space, starting at
! index is in s. If ss is found, is is set to (index of last character of ss in
! s) + 1; else is is set to 0. If is is out of range on input, routine
! terminates with is = -1.
character(*), intent(in) :: s   ! input string
integer, intent(inout) :: is    ! on input: starting index for search for ss in
                                ! s on output: (index of last character of ss in
                                ! s) + 1
character(*), intent(out) :: ss ! first substring in s, starting from index is
character(len(s)+1) :: t        ! temporary string to facilitate search
integer i, i1, i2
logical prevwhite, curwhite
if (is <= 0 .or. is > len(s)) then
    ss = ""; is = -1; return
end if
t = s // " "
if (is == 1) then
    prevwhite = .true.
else
    prevwhite = whitechar(t(is-1:is-1))
end if
i1 = 0; i2 = 0
do i = is, len(t)
    curwhite = whitechar(t(i:i))
    if (prevwhite .and. .not. curwhite) i1 = i   ! beginning of substring
    if (i1>0 .and. curwhite) then                ! end of substring
        i2 = i-1; exit
    end if
    prevwhite=curwhite
end do
if (i2 > 0) then
    ss = t(i1:i2); is = i2+1
else
    ss = ""; is = 0
end if
end subroutine

integer function newunit(unit) result(n)
! Returns lowest i/o unit number not in use (to be used in older compilers).
!
! Starting at 10 to avoid lower numbers which are sometimes reserved.
! Note: largest valid unit number may be system-dependent.
!
! Arguments
! ---------
!
! If present, the new unit will be returned into it
integer, intent(out), optional :: unit
!
! Example
! -------
!
! integer :: u
! open(newunit(u), file="log.txt", status="old")
! read(u, *) a, b
! close(u)
!
! In new compilers, just use the "newunit" keyword argument:
!
! integer :: u
! open(newunit=u, file="log.txt", status="old")
! read(u, *) a, b
! close(u)

logical inuse
integer, parameter :: nmin=10   ! avoid lower numbers which are sometimes reserved
integer, parameter :: nmax=999  ! may be system-dependent
do n = nmin, nmax
    inquire(unit=n, opened=inuse)
    if (.not. inuse) then
        if (present(unit)) unit=n
        return
    end if
end do
call stop_error("newunit ERROR: available unit not found.")
end function

subroutine stop_error(msg)
! Aborts the program with nonzero exit code
!
! The statement "stop msg" will return 0 exit code when compiled using
! gfortran. stop_error() uses the statement "stop 1" which returns an exit code
! 1 and a print statement to print the message.
!
! Example
! -------
!
! call stop_error("Invalid argument")

character(len=*) :: msg ! Message to print on stdout
print *, msg
stop 1
end subroutine

subroutine loadtxt(filename, d)
! Loads a 2D array from a text file.
!
! Arguments
! ---------
!
! Filename to load the array from
character(len=*), intent(in) :: filename
! The array 'd' will be automatically allocated with the correct dimensions
real(dp), allocatable, intent(out) :: d(:, :)
!
! Example
! -------
!
! real(dp), allocatable :: data(:, :)
! call loadtxt("log.txt", data)  ! 'data' will be automatically allocated
!
! Where 'log.txt' contains for example::
!
!     1 2 3
!     2 4 6
!     8 9 10
!     11 12 13
!     ...
!
character :: c
integer :: s, ncol, nrow, ios, i
logical :: lastwhite
real(dp) :: r

open(newunit(s), file=filename, status="old")

! determine number of columns
ncol = 0
lastwhite = .true.
do
   read(s, '(a)', advance='no', iostat=ios) c
   if (ios /= 0) exit
   if (lastwhite .and. .not. whitechar(c)) ncol = ncol + 1
   lastwhite = whitechar(c)
end do

rewind(s)

! determine number or rows
nrow = 0
do
   read(s, *, iostat=ios) r
   if (ios /= 0) exit
   nrow = nrow + 1
end do

rewind(s)

allocate(d(nrow, ncol))
do i = 1, nrow
    read(s, *) d(i, :)
end do
close(s)
end subroutine

subroutine savetxt(filename, d)
! Saves a 2D array into a textfile.
!
! Arguments
! ---------
!
character(len=*), intent(in) :: filename  ! File to save the array to
real(dp), intent(in) :: d(:, :)           ! The 2D array to save
!
! Example
! -------
!
! real(dp) :: data(3, 2)
! call savetxt("log.txt", data)

integer :: s, i
open(newunit(s), file=filename, status="replace")
do i = 1, size(d, 1)
    write(s, *) d(i, :)
end do
close(s)
end subroutine

subroutine arange(a, b, dx, u)
! Returns an array u = [a, a+dx, a+2*dx, ..., b-dx]
!
! Arguments
! ---------
!
real(dp), intent(in) :: a, b, dx
real(dp), allocatable, intent(out) :: u(:)
!
! Example
! -------
!
! real(dp), allocatable :: u(:)
! call arange(1, 5, 1, u)   ! u = [1, 2, 3, 4]
integer :: n, i
n = int((b-a) / dx)
allocate(u(n))
do i = 1, n
    u(i) = a + (i-1)*dx
end do
end subroutine

function zeros(n) result(x)
integer, intent(in) :: n
real(dp) :: x(n)
x = 0
end function

subroutine assert(condition)
! If condition == .false., it aborts the program.
!
! Arguments
! ---------
!
logical, intent(in) :: condition
!
! Example
! -------
!
! call assert(a == 5)

if (.not. condition) call stop_error("Assert failed.")
end subroutine

pure integer function str_int_len(i) result(sz)
! Returns the length of the string representation of 'i'
integer, intent(in) :: i
integer, parameter :: MAX_STR = 100
character(MAX_STR) :: s
! If 's' is too short (MAX_STR too small), Fortan will abort with:
! "Fortran runtime error: End of record"
write(s, '(i0)') i
sz = len_trim(s)
end function

pure function str_int(i) result(s)
! Converts integer "i" to string
integer, intent(in) :: i
character(len=str_int_len(i)) :: s
write(s, '(i0)') i
end function

pure integer function str_real_len(r, fmt) result(sz)
! Returns the length of the string representation of 'i'
real(dp), intent(in) :: r
character(len=*), intent(in) :: fmt
integer, parameter :: MAX_STR = 100
character(MAX_STR) :: s
! If 's' is too short (MAX_STR too small), Fortan will abort with:
! "Fortran runtime error: End of record"
write(s, fmt) r
sz = len_trim(s)
end function

pure function str_real(r) result(s)
! Converts the real number "r" to string with 7 decimal digits.
real(dp), intent(in) :: r
character(len=*), parameter :: fmt="(f0.6)"
character(len=str_real_len(r, fmt)) :: s
write(s, fmt) r
end function

pure function str_real_n(r, n) result(s)
! Converts the real number "r" to string with 'n' decimal digits.
real(dp), intent(in) :: r
integer, intent(in) :: n
character(len=str_real_len(r, "(f0." // str_int(n) // ")")) :: s
write(s, "(f0." // str_int(n) // ")") r
end function

pure integer function strfmt_int_len(fm, i) result(sz)
! Returns the length of the string representation of 'i'
character(len=*), intent(in) :: fm
integer, intent(in) :: i
integer, parameter :: MAX_STR = 100
character(MAX_STR) :: s
! If 's' is too short (MAX_STR too small), Fortan will abort with:
! "Fortran runtime error: End of record"
write(s, fm) i
sz = len_trim(s)
end function

pure function strfmt_int(fm, i) result(s)
! Converts integer "i" to string
character(len=*), intent(in) :: fm
integer, intent(in) :: i
character(len=strfmt_int_len(fm, i)) :: s
write(s, fm) i
end function

pure integer function strfmt_real_len(fm, r) result(sz)
! Returns the length of the string representation of 'r'
character(len=*), intent(in) :: fm
real(dp), intent(in) :: r
integer, parameter :: MAX_STR = 100
character(MAX_STR) :: s
! If 's' is too short (MAX_STR too small), Fortan will abort with:
! "Fortran runtime error: End of record"
write(s, fm) r
sz = len_trim(s)
end function

pure function strfmt_real(fm, r) result(s)
! Converts real "r" to string
character(len=*), intent(in) :: fm
real(dp), intent(in) :: r
character(len=strfmt_real_len(fm, r)) :: s
write(s, fm) r
end function

subroutine init_random()
! Initializes the random number generator based on the system's time.
integer :: i, n, clock
integer, allocatable :: seed(:)
call random_seed(size=n)
allocate(seed(n))
call system_clock(count=clock)
seed = clock + 37 * [(i - 1, i = 1, n)]
call random_seed(put=seed)
end subroutine

function linspace(a, b, n) result(s)
real(dp), intent(in) :: a, b
integer, intent(in) :: n
real(dp) :: s(n)
s = mesh_exp(a, b, 1.0_dp, n-1)
end function

function mesh_exp(r_min, r_max, a, N) result(mesh)
! Generates exponential mesh of N elements on [r_min, r_max]
!
! Arguments
! ---------
!
! The domain [r_min, r_max], the mesh will contain both endpoints:
real(dp), intent(in) :: r_min, r_max
!
! The fraction of the rightmost vs. leftmost elements of the mesh (for a > 1
! this means the "largest/smallest"); The only requirement is a > 0. For a == 1
! a uniform mesh will be returned:
real(dp), intent(in) :: a
!
! The number of elements in the mesh:
integer, intent(in) :: N
!
! Returns
! -------
!
! The generated mesh:
real(dp) :: mesh(N+1)
!
! Note: Every exponential mesh is fully determined by the set of parameters
! (r_min, r_max, a, N). Use the get_mesh_exp_params() subroutine to obtain them
! from the given mesh.
!
! Example
! -------
!
! real(dp) :: r(11)
! r = mesh_exp(0._dp, 50._dp, 1e9_dp, 10)

integer :: i
real(dp) :: alpha, beta
if (a < 0) then
    call stop_error("mesh_exp: a > 0 required")
else if (abs(a - 1) < 1e-16_dp) then
    alpha = (r_max - r_min) / N
    do i = 1, N+1
        mesh(i) = alpha * (i-1.0_dp) + r_min
    end do
else
    if (N > 1) then
        beta = log(a) / (N-1)
        alpha = (r_max - r_min) / (exp(beta*N) - 1)
        do i = 1, N+1
            mesh(i) = alpha * (exp(beta*(i-1)) - 1) + r_min
        end do
    else if (N == 1) then
        mesh(1) = r_min
        mesh(2) = r_max
    else
        call stop_error("mesh_exp: N >= 1 required")
    end if
end if
end function

real(dp) function clock() result(r)
use openmp, only: omp_get_wtime, with_openmp
if (with_openmp()) then
    r = omp_get_wtime()
else
    call cpu_time(r)
end if
end function


integer function get_int_arg(i) result(r)
integer, intent(in) :: i
integer, parameter :: maxlen=32
character(len=maxlen) :: arg
integer :: s
if (.not. (0 < i .and. i <= command_argument_count())) then
    call stop_error("get_int_arg: `i` must satisfy `0 < i <= arg_count`.")
end if
call get_command_argument(i, arg, status=s)
if (s == -1) then
    call stop_error("get_int_arg: Argument too long, increase `maxlen`.")
else if (s > 0) then
    call stop_error("get_int_arg: Argument retrieval failed.")
else if (s /= 0) then
    call stop_error("get_int_arg: Unknown error.")
end if
read(arg, *, iostat=s) r
if (s /= 0) then
    call stop_error("get_int_arg: Failed to convert an argument to integer.")
end if
end function

real(dp) function get_float_arg(i) result(r)
integer, intent(in) :: i
integer, parameter :: maxlen=32
character(len=maxlen) :: arg
integer :: s
if (.not. (0 < i .and. i <= command_argument_count())) then
    call stop_error("get_int_arg: `i` must satisfy `0 < i <= arg_count`.")
end if
call get_command_argument(i, arg, status=s)
if (s == -1) then
    call stop_error("get_int_arg: Argument too long, increase `maxlen`.")
else if (s > 0) then
    call stop_error("get_int_arg: Argument retrieval failed.")
else if (s /= 0) then
    call stop_error("get_int_arg: Unknown error.")
end if
read(arg, *, iostat=s) r
if (s /= 0) then
    call stop_error("get_int_arg: Failed to convert an argument to a float.")
end if
end function

subroutine allocate_mold_real3d(A, B)
real(dp), allocatable, intent(out) :: A(:,:,:)
real(dp), intent(in) :: B(:,:,:)
! Equivalent to allocate(A, mold=B). Use this with compilers that do not
! support the F2008 syntax yet.
allocate(A(size(B,1), size(B,2), size(B,3)))
end subroutine

subroutine allocate_mold_complex3d(A, B)
complex(dp), allocatable, intent(out) :: A(:,:,:)
complex(dp), intent(in) :: B(:,:,:)
! Equivalent to allocate(A, mold=B). Use this with compilers that do not
! support the F2008 syntax yet.
allocate(A(size(B,1), size(B,2), size(B,3)))
end subroutine

end module
