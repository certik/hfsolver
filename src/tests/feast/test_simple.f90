program test_simple
use types, only: dp
use utils, only: stop_error, assert
use linalg, only: eye
use linalg_feast, only: eigh
implicit none

real(dp), parameter :: eps = 1e-9_dp
real(dp) :: A(2, 2) = reshape([2, -1, -1, 2], [2, 2])
real(dp) :: Emin=-10, Emax=10
integer :: M0=2 ! (Initial) subspace dimension
real(dp), allocatable :: lam(:), c(:, :), r(:)
real(dp) :: norm
integer :: i

allocate(lam(M0))
allocate(r(M0))
allocate(c(2,M0))

call eigh(A, Emin, Emax, M0, lam, c)

print *, 'Eigenvalues'
do i = 1, M0
    print *, i, lam(i)
end do
print *,'Eigenvectors'
do i=1, M0
    print *,i,"(",c(1,i),c(2,i),")"
end do

call assert(all(abs(lam - [1, 3]) < eps))
do i = 1, 2
    ! Test that c(:, i) are eigenvectors:
    r = matmul(A-lam(i)*eye(2), c(:, i))
    call assert(sqrt(dot_product(r, r)) < eps)
    ! Test that eigenvectors are properly normalized:
    norm = dot_product(c(:, i), c(:, i))
    call assert(abs(norm - 1) < eps)
end do

end program
