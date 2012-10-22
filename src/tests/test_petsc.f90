program test_petsc
use types, only: dp
use constants, only: i_
use petsc_, only: petsc_init, petsc_finalize, solve
implicit none

integer, parameter :: n = 10
complex(dp) :: A(n, n), b(n), x(n)
integer :: i, ierr

A = 0
A(1, 1:2) = [2, -1]
forall(i = 2:n-1) A(i, i-1:i+1) = [-1, 2, -1]
A(n, n-1:n) = [-1, 2]

x = 1+2*i_
b = matmul(A, x)

ierr = petsc_init()
ierr = solve(n, A, b, x)
print *, x
ierr = petsc_finalize()
end
