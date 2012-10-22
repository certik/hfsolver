module sparse

! These routines are following similar API and algorithms as routines in
! scipy.sparse.sparsetools

use types, only: dp
use sorting, only: argsort
implicit none
private
public coo2dense, dense2coo, getnnz, coo2csr, coo2csc, &
    csr_has_canonical_format, csr_sum_duplicates, csr_sort_indices, &
    coo2csr_canonical, csr_matvec

contains

! Dense

subroutine dense2coo(B, Ai, Aj, Ax)
real(dp), intent(in) :: B(:, :)
integer, intent(out) :: Ai(:), Aj(:)
real(dp), intent(out) :: Ax(:)
integer :: i, j, idx
idx = 1
do j = 1, size(B, 2)
    do i = 1, size(B, 1)
        if (B(i, j) == 0) cycle
        Ai(idx) = i
        Aj(idx) = j
        Ax(idx) = B(i, j)
        idx = idx + 1
    end do
end do
end subroutine

integer function getnnz(B) result(nnz)
real(dp), intent(in) :: B(:, :)
integer :: i, j
nnz = 0
do j = 1, size(B, 2)
    do i = 1, size(B, 1)
        if (B(i, j) == 0) cycle
        nnz = nnz + 1
    end do
end do
end function

! COO

subroutine coo2dense(Ai, Aj, Ax, B)
integer, intent(in) :: Ai(:), Aj(:)
real(dp), intent(in) :: Ax(:)
real(dp), intent(out) :: B(:, :)
integer :: n
B = 0
do n = 1, size(Ai)
    B(Ai(n), Aj(n)) = B(Ai(n), Aj(n)) + Ax(n)
end do
end subroutine

subroutine coo2csr(Ai, Aj, Ax, Bp, Bj, Bx)
! Converts from COO (Ai, Aj, Ax) into CSR (Bp, Bj, Bx)
! Row and column indices are *not* assumed to be ordered.
! Duplicate entries are carried over to the CSR representation.
integer, intent(in) :: Ai(:), Aj(:)
real(dp), intent(in) :: Ax(:)
integer, intent(out) :: Bp(:), Bj(:)
real(dp), intent(out) :: Bx(:)
integer :: n, i, n_row, nnz, cumsum, temp, row, dest
n_row = size(Bp)-1
nnz = size(Ai)
Bp = 0
do n = 1, nnz
    Bp(Ai(n)) = Bp(Ai(n)) + 1
end do
cumsum = 1
do i = 1, n_row
    temp = Bp(i)
    Bp(i) = cumsum
    cumsum = cumsum + temp
end do
do n = 1, nnz
    row = Ai(n)
    dest = Bp(row)
    Bj(dest) = Aj(n)
    Bx(dest) = Ax(n)
    Bp(row) = Bp(row) + 1
end do
Bp(2:) = Bp(:n_row)
Bp(1) = 1
end subroutine

subroutine coo2csc(Ai, Aj, Ax, Bp, Bi, Bx)
! Converts from COO (Ai, Aj, Ax) into CSC (Bp, Bi, Bx)
! Row and column indices are *not* assumed to be ordered.
! Duplicate entries are carried over to the CSC representation.
integer, intent(in) :: Ai(:), Aj(:)
real(dp), intent(in) :: Ax(:)
integer, intent(out) :: Bp(:), Bi(:)
real(dp), intent(out) :: Bx(:)
! Calculate CSR of the transposed matrix:
call coo2csr(Aj, Ai, Ax, Bp, Bi, Bx)
end subroutine

subroutine coo2csr_canonical(Ai, Aj, Ax, Bp, Bj, Bx)
! Converts from COO (Ai, Aj, Ax) into CSR (Bp, Bj, Bx)
! Row and column indices are *not* assumed to be ordered.
! Duplicate entries are summed up and the indices are ordered.
integer, intent(in) :: Ai(:), Aj(:)
real(dp), intent(in) :: Ax(:)
integer, allocatable, intent(out) :: Bp(:), Bj(:)
real(dp), allocatable, intent(out) :: Bx(:)
integer, allocatable :: Bj_(:)
real(dp), allocatable :: Bx_(:)
integer :: nnz
allocate(Bp(maxval(Ai)+1), Bj_(size(Ai)), Bx_(size(Ai)))
call coo2csr(Ai, Aj, Ax, Bp, Bj_, Bx_)
call csr_sort_indices(Bp, Bj_, Bx_)
call csr_sum_duplicates(Bp, Bj_, Bx_)
nnz = Bp(size(Bp))-1
allocate(Bj(nnz), Bx(nnz))
Bj = Bj_(:nnz)
Bx = Bx_(:nnz)
end subroutine

! CSR

logical function csr_has_canonical_format(Ap, Aj) result(r)
! Determine whether the matrix structure is canonical CSR.
! Canonical CSR implies that column indices within each row
! are (1) sorted and (2) unique.  Matrices that meet these
! conditions facilitate faster matrix computations.
integer, intent(in) :: Ap(:), Aj(:)
integer :: i, j
r = .false.
do i = 1, size(Ap)-1
    if (Ap(i) > Ap(i+1)) return
    do j = Ap(i)+1, Ap(i+1)-1
        if (Aj(j-1) >= Aj(j)) return
    end do
end do
r = .true.
end function

subroutine csr_sum_duplicates(Ap, Aj, Ax)
! Sum together duplicate column entries in each row of CSR matrix A
! The column indicies within each row must be in sorted order.
! Explicit zeros are retained.
! Ap, Aj, and Ax will be modified *inplace*
integer, intent(inout) :: Ap(:), Aj(:)
real(dp), intent(inout) :: Ax(:)
integer :: nnz, r1, r2, i, j, jj
real(dp) :: x
nnz = 1
r2 = 1
do i = 1, size(Ap) - 1
    r1 = r2
    r2 = Ap(i+1)
    jj = r1
    do while (jj < r2)
        j = Aj(jj)
        x = Ax(jj)
        jj = jj + 1
        do while (jj < r2)
            if (Aj(jj) == j) then
                x = x + Ax(jj)
                jj = jj + 1
            else
                exit
            end if
        end do
        Aj(nnz) = j
        Ax(nnz) = x
        nnz = nnz + 1
    end do
    Ap(i+1) = nnz
end do
end subroutine

subroutine csr_sort_indices(Ap, Aj, Ax)
! Sort CSR column indices inplace
integer, intent(inout) :: Ap(:), Aj(:)
real(dp), intent(inout) :: Ax(:)
integer :: i, r1, r2, l, idx(size(Aj))
do i = 1, size(Ap)-1
    r1 = Ap(i)
    r2 = Ap(i+1)-1
    l = r2-r1+1
    idx(:l) = argsort(Aj(r1:r2))
    Aj(r1:r2) = Aj(r1+idx(:l)-1)
    Ax(r1:r2) = Ax(r1+idx(:l)-1)
end do
end subroutine

subroutine csr_matvec(Ap, Aj, Ax, x, y)
! Compute y = A*x for CSR matrix A and dense vectors x, y
integer, intent(in) :: Ap(:), Aj(:)
real(dp), intent(in) :: Ax(:), x(:)
real(dp), intent(out) :: y(:)
integer :: i
forall(i=1:size(Ap)-1) y(i) = sum(Ax(Ap(i):Ap(i+1)-1) * x(Aj(Ap(i):Ap(i+1)-1)))
end subroutine

end module
