program test_fourier
use types, only: dp
use fourier, only: dft, idft
use utils, only: assert
use constants, only: i_
implicit none

real(dp), allocatable :: x(:)
complex(dp), allocatable :: xdft(:)
integer :: n, i

call assert(all(abs(dft([1._dp]) - [1]) < 1e-12_dp))
call assert(all(abs(idft([1._dp+0*i_]) - [1]) < 1e-12_dp))
call assert(all(abs(dft([1._dp, 2._dp]) - [3, -1]) < 1e-12_dp))
call assert(all(abs(idft([3._dp+0*i_, -1._dp + 0*i_]) - [1, 2]) < 1e-12_dp))

allocate(x(4), xdft(4))
x = [1, 2, 3, 4]
xdft = [10 + 0*i_, -2+2*i_, -2+0*i_, -2-2*i_]
call assert(all(abs(dft(x) - xdft) < 1e-12_dp))
call assert(all(abs(idft(xdft) - x) < 1e-12_dp))
deallocate(x, xdft)

allocate(x(5))
x = [1, 2, 3, 4, 5]
call assert(all(abs(idft(dft(x)) - x) < 1e-12_dp))
deallocate(x)

n = 400
allocate(x(n))
forall(i = 1:n) x(i) = i
call assert(all(abs(idft(dft(x)) - x) < 1e-10_dp))
deallocate(x)

end program
