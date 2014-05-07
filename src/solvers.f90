module solvers
use types, only: dp
use lapack, only: dlamch, dsygvx, dggev, dsysv, dsygvd
use utils, only: stop_error
use sorting, only: sortpairs
implicit none
private
public solve_eig_range, spurious, filter_spurious_states, count_nodes, &
    solve_eig_range2, solve_sym, solve_eig_range3

CONTAINS

subroutine solve_eig_range(Am, Bm, Emin, Emax, lam, c)
! solves system equations
! simple dense solver for now, sparse solver better suited for low-order FE/SE ...
real(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam Bm c
real(dp), intent(in) :: Bm(:,:)   ! RHS matrix: Am c = lam Bm c
real(dp), intent(in) :: Emin, Emax
real(dp), allocatable, intent(out) :: lam(:)   ! eigenvalues: Am c = lam Bm c
real(dp), allocatable, intent(out) :: c(:,:)   ! eigenvectors: Am c = lam Bm c; c(i,j) = ith component of jth vec.
integer i, n
! lapack variables
integer lwork, info, m
integer, allocatable:: iwork(:), ifail(:)
real(dp) abstol
real(dp), allocatable:: Amt(:,:), Bmt(:,:), lamt(:), work(:)
real(dp), allocatable:: z(:,:)

! solve
n=size(Am,1)
lwork=8*n
allocate(Amt(n,n),Bmt(n,n),lamt(n),work(lwork),iwork(5*n),ifail(n),z(n,n))
Amt=Am; Bmt=Bm  ! Amt,Bmt temporaries overwritten by dsygvx
abstol=2*dlamch('S')
call dsygvx(1,'V','V','L',n,Amt,n,Bmt,n,Emin,Emax,0,0,abstol,m, &
    lamt,z,n,work,lwork,iwork,ifail,info)
if (info/=0) then
    print *, "dsygvx returned info =", info
    if (info <= n) then
        print *, info, "eigenvectors failed to converge. Their indices are:"
        print *, ifail
    else
        print *, "The leading minor of order ", info-n, &
            "of B is not positive definite. The factorization of B could ", &
            "not be completed and no eigenvalues or eigenvectors were computed."
    end if
    call stop_error('DSYGVX ERROR')
end if
! Copy the found eigenvalues/vectors into lam/c
allocate(lam(m), c(n, m))
lam(:) = lamt(:m)
c(:, :) = z(:, :m)
do i = 1, m ! normalize such that initial values positive
   if (c(2,i)<0) c(:,i)=-c(:,i)
end do
!write(*,*) "Eigenvector residuals ||Am c - lam Bm c||: "
!allocate(r(n))
!do i = 1, m
!   r=matmul(Am-lam(i)*Bm,c(:,i))
!   write(*,'(1x,i5,a,es18.11)') i, ": ", sqrt(dot_product(r,r))
!   print *, dot_product(c(:, i), matmul(Bm, c(:, i)))
!end do
end subroutine

subroutine solve_eig_range2(Am, Bm, Emin, Emax, lam, c)
! solves generalized nonsymmetric eigenvalue problem
real(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam Bm c
real(dp), intent(in) :: Bm(:,:)   ! RHS matrix: Am c = lam Bm c
real(dp), intent(in) :: Emin, Emax
real(dp), allocatable, intent(out) :: lam(:)   ! eigenvalues: Am c = lam Bm c
real(dp), allocatable, intent(out) :: c(:,:)   ! eigenvectors: Am c = lam Bm c; c(i,j) = ith component of jth vec.
integer i, j, n
! lapack variables
integer lwork, info
real(dp), allocatable:: Amt(:,:), Bmt(:,:), work(:)
real(dp), allocatable:: r(:)
real(dp), allocatable :: vl(:, :)
real(dp), allocatable :: vr(:, :), alphar(:), alphai(:), beta(:)
real(dp), allocatable :: lambda(:)
integer, allocatable :: sel(:)
!complex(dp), parameter :: J = (0, 1)

! solve
n=size(Am,1)
lwork=8*n
allocate(Amt(n,n), Bmt(n,n), work(lwork), vr(n, n), alphar(n), alphai(n), &
    beta(n), vl(n, n))
Amt=Am; Bmt=Bm  ! Amt,Bmt temporaries overwritten by dsygvx
call dggev('N', 'V', n, Amt, n, Bmt, n, alphar, alphai, beta, vl, n, &
    vr, n, work, lwork, info)
if (info/=0) then
    print *, "dsygvx returned info =", info
    if (info <= n) then
        print *, "The", info, " iteration failed."
    else
        print *, "Other error"
    end if
    call stop_error('DGGEV ERROR')
end if
allocate(lambda(n), sel(n))
where (abs(beta) > 1e-16_dp) lambda = alphar / beta
j = 0
do i = 1, n
    sel(i) = 0
    if (alphai(i) == 0) then
        if (beta(i) == 0) then
            ! Maybe we can just ignore the beta==0 cases
            call stop_error("Got beta==0")
        end if
        if (lambda(i) > Emin .and. lambda(i) < Emax) then
            j = j + 1
            sel(i) = j
        end if
    end if
end do
allocate(lam(maxval(sel)), c(n, maxval(sel)))
do i = 1, size(lambda)
    if (sel(i) /= 0) then
        lam(sel(i)) = lambda(i)
        c(:, sel(i)) = vr(:, i)
    end if
end do
call sortpairs(lam, c)
do i = 1, size(lam) ! normalize such that initial values positive
   if (c(2,i)<0) c(:,i)=-c(:,i)
end do
! Normalize the eigenvectors:
do i = 1, size(lam)
   c(:, i) = c(:, i) / sqrt(dot_product(c(:, i), matmul(Bm, c(:, i))))
end do
write(*,*) "Eigenvector residuals ||Am c - lam Bm c||: "
allocate(r(n))
do i = 1, size(lam)
   r=matmul(Am-lam(i)*Bm,c(:,i))
   write(*,'(1x,i5,a,es18.11)') i, ": ", sqrt(dot_product(r,r))
   !print *, dot_product(c(:, i), matmul(Bm, c(:, i)))
end do
end subroutine

subroutine solve_eig_range3(Am, Bm, Emin, Emax, lam, c)
! solves system equations
! simple dense solver for now, sparse solver better suited for low-order FE/SE ...
real(dp), intent(in) :: Am(:,:)   ! LHS matrix: Am c = lam Bm c
real(dp), intent(in) :: Bm(:,:)   ! RHS matrix: Am c = lam Bm c
real(dp), intent(in) :: Emin, Emax
real(dp), allocatable, intent(out) :: lam(:)   ! eigenvalues: Am c = lam Bm c
real(dp), allocatable, intent(out) :: c(:,:)   ! eigenvectors: Am c = lam Bm c; c(i,j) = ith component of jth vec.
integer i, n
! lapack variables
integer lwork, liwork, info, m, m1, m2
integer, allocatable:: iwork(:)
real(dp), allocatable:: Amt(:,:), Bmt(:,:), lamt(:), work(:)
real(dp), allocatable:: r(:)

! solve
n=size(Am,1)
lwork=1+6*n+2*n**2
liwork=3+5*n
allocate(Amt(n,n),Bmt(n,n),lamt(n),work(lwork),iwork(liwork))
Amt=Am; Bmt=Bm  ! Amt,Bmt temporaries overwritten by dsygvx
call dsygvd(1,'V','L',n,Amt,n,Bmt,n,lamt,work,lwork,iwork,liwork,info)
if (info/=0) then
    print *, "dsygvd returned INFO =", info
    if (info < 0) then
        print *, "the", -info, "-th argument had an illegal value"
    else if (info <= n) then
        print *, "the algorithm failed to compute an eigenvalue while working"
        print *, "on the submatrix lying in rows and columns", 1.0_dp*info/(n+1)
        print *, "through", mod(info, n+1)
    else
        print *, "The leading minor of order ", info-n, &
            "of B is not positive definite. The factorization of B could ", &
            "not be completed and no eigenvalues or eigenvectors were computed."
    end if
    call stop_error('DSYGVD ERROR')
end if
! TODO: make sure things work with regards m1==n and m2==n:
m1 = 1
do while (lamt(m1) < Emin .and. m1 < n)
    m1 = m1 + 1
end do
! m1 is the first eigenvalue in the range
m2 = m1
do while (lamt(m2) <= Emax .and. m2 < n)
    m2 = m2 + 1
end do
m2 = m2 - 1
! m2 is the last eigenvalue in the range
m = m2-m1+1 ! number of eigenvalues in the range
! Copy the found eigenvalues/vectors into lam/c
allocate(lam(m), c(n, m))
lam(:) = lamt(m1:m2)
c(:, :) = Amt(:, m1:m2)
do i = 1, m ! normalize such that initial values positive
   if (c(2,i)<0) c(:,i)=-c(:,i)
end do
write(*,*) "Eigenvector residuals ||Am c - lam Bm c||: "
allocate(r(n))
do i = 1, m
   r=matmul(Am-lam(i)*Bm,c(:,i))
   write(*,'(1x,i5,a,es18.11)') i, ": ", sqrt(dot_product(r,r))
   print *, dot_product(c(:, i), matmul(Bm, c(:, i)))
end do
end subroutine

integer function count_nodes(u) result(nn)
! Returns the number of nodes of the function
real(dp), intent(in) :: u(:) ! State "u"

integer :: i
logical :: s, last_s, ok
real(dp) :: thres
thres = maxval(abs(u)) * 1e-5_dp
nn = 0
ok = .false.
do i = 1, size(u)
    ! skip small "tails" of the wavefunction whose oscillations we don't want to
    ! count
    if (abs(u(i)) < thres) cycle
    s = u(i) > 0
    if (ok) then
        if (s .neqv. last_s) nn = nn + 1
    end if
    last_s = s
    ok = .true.
end do
end function

logical function spurious_old2(u)
! Returns .true. if the state "u" is spurious, otherwise .false.
! The state is declared spurious if it has (on average) less nodes per
! oscillation than allowed (see the nn_per_oscillation below), determined by
! counting the number of oscillations and comparing to the length of the vector.
real(dp), intent(in) :: u(:) ! State "u"
! The minimum allowed nodes per oscillation
integer, parameter :: nn_per_oscillation = 4

real(dp) :: osc(size(u)-1)
integer :: i, nn
logical :: s, last_s
real(dp) :: thres
integer :: last_oscillation
osc(:) = u(2:) - u(:size(u)-1) ! Calculate differences
thres = (maxval(u) - minval(u)) * 1e-5_dp ! Determine treshold
!print *, "threshold:", thres
nn = 0
last_s = osc(1) > 0
!print *, "u:"
!print *, u
!print *, "osc:"
!print *, osc
last_oscillation = 1
do i = 1, size(osc)
    ! skip small "tails" of the wavefunction whose oscillations we don't want to
    ! count
    if (abs(osc(i)) < thres) cycle
    s = osc(i) > 0
    if (s .neqv. last_s) then
        nn = nn + 1
        !if (i==99) then
        !    print *, u(i-3:)
        !end if
!        print *, "found oscillation at:", i, size(u), &
!            abs(u(i) - u(last_oscillation))
        if (i - last_oscillation < nn_per_oscillation) then
            spurious_old2 = .true.
            return
        end if
        last_oscillation = i
    end if
    last_s = s
end do
if (nn > size(u) / nn_per_oscillation) then
    spurious_old2 = .true.
else
    spurious_old2 = .false.
end if
!print *, "oscillations:", nn, size(u) / nn_per_oscillation
!print *, maxval(abs(u)), u(:10)
!print *, thres, osc(:10)
!spurious = .false.
end function

logical function spurious(u, ib)
! Returns .true. if the state "u" is spurious, otherwise .false.
! The state is declared spurious if it has (on average) less nodes per
! oscillation than allowed (see the nn_per_oscillation below), determined by
! counting the number of oscillations and comparing to the length of the vector.
real(dp), intent(in) :: u(:) ! State "u"
integer, intent(in) :: ib(:, :)   ! basis connectivity: ib(i,j) = index of basis
   ! function associated with local basis function i of element j. 0 = no
   ! associated basis fn.
! The minimum allowed nodes per oscillation
!integer, parameter :: nn_per_oscillation = 4

real(dp) :: osc(size(u) - 1)
integer :: i, ie, ip, nn
logical :: s, last_s
real(dp) :: thres
integer :: last_oscillation
osc(:) = u(2:) - u(:size(u)-1) ! Calculate differences
thres = (maxval(u) - minval(u)) * 1e-5_dp ! Determine treshold
!print *, "threshold:", thres
nn = 0
last_s = osc(1) > 0
!print *, "u:"
!print *, u
!print *, "osc:"
!print *, osc
if (size(ib, 1)-1 == 1) call stop_error("spurious check doesn't work for p=1")
last_oscillation = 1
do ie = 1, size(ib, 2)
    do ip = 1, size(ib, 1)
        i = ib(ip, ie)
        if (i == 0) cycle
        if (i < 0) cycle
        if (i > size(osc)) exit
        s = osc(i) > 0
        if (s .neqv. last_s) then
            if (abs(u(i)-u(last_oscillation)) < thres) then
                last_oscillation = i
                last_s = s
                cycle
            end if
            nn = nn + 1
            last_oscillation = i
        end if
        last_s = s
    end do
    !print *, "nn = ", nn, size(ib, 1), (size(ib, 1)-1)/3 + 1
    if (nn > (size(ib, 1)-1)/3 + 1) then
        print *, "spurious, nodes:", nn, "element:", ie
        spurious = .true.
        return
    end if
    nn = 0
end do
spurious = .false.
end function

subroutine filter_spurious_states(ib, lam, u, lam2, u2, skipped)
integer, intent(in) :: ib(:,:,:)   ! basis connectivity: ib(i,j) = index of basis
   ! function associated with local basis function i of element j. 0 = no
   ! associated basis fn.
real(dp), intent(in) :: u(:, :), lam(:)
real(dp), intent(out), allocatable :: u2(:, :), lam2(:)
logical, intent(out) :: skipped
integer :: m(size(lam))
integer :: i, j, c, c_begin, c_end, n_eig
logical :: sp
j = 0
print *, "Filtering high oscillatory states:"
do i = 1, size(lam)
    ! Loop over the vector components:
    sp = .false.
    c_begin = 1
    do c = 1, size(ib, 3)
        c_end = maxval(ib(:, :, c))
        if (spurious(u(c_begin:c_end, i), ib(:, :, c))) then
            sp = .true.
            exit
        end if
        c_begin = c_end+1
    end do
    if (sp) then
        m(i) = 0
    else
        j = j + 1
        m(i) = j
    end if
    print "(i4,i4,f12.5)", i, m(i), lam(i)
end do
n_eig = maxval(m)
allocate(u2(size(u, 1), n_eig), lam2(n_eig))
do i = 1, size(lam)
    if (m(i) /= 0) then
        u2(:, m(i)) = u(:, i)
        lam2(m(i)) = lam(i)
    end if
end do
! If we skipped any state (we filtered some in between), then these two numbers
! will mismatch:
!skipped = m(n_eig) /= n_eig
! But we actually want to see, whether we filtered anything at all:
if (size(m) == 0) then
    skipped = .false.
else
    skipped = m(size(m)) /= size(m)
end if
end subroutine

function solve_sym(Am, bv) result(c)
! solves symmetric dense system of equations
real(dp), intent(in) :: Am(:,:)   ! system matrix: Am c = bv
real(dp), intent(in) :: bv(:)     ! source vector: Am c = bv
real(dp) :: c(size(bv))     ! solution vector: Am c = bv
!real(dp) :: r(size(bv))
integer :: n
! lapack variables
integer :: lwork, info
integer, allocatable :: ipiv(:)
real(dp), allocatable :: Amt(:,:),bm(:,:),work(:)

n = size(c)
lwork = n
allocate(Amt(n,n), bm(n,1), ipiv(n), work(lwork))
Amt=Am; bm(:,1)=bv   ! temporaries for dsysv
call dsysv('L', n, 1, Amt, n, ipiv, bm, n, work, lwork, info)
if (info < 0) then
    print *, "The", -info, "-th argument had illegal value"
    call stop_error('DSYSV ERROR.')
end if
if (info > 0) then
    print *, "D(", info, ",", info, ") is exactly zero."
    print *, "The factorization has been completed, but the block diagonal"
    print *, "matrix D is exactly singular, so the solution could not be"
    print *, "computed."
    call stop_error('DSYSV ERROR.')
end if
c = bm(:, 1)
! error
!r=matmul(Am, c)-bv
!write(*,'(1x,a,es18.11)') "Solution vector residual ||Am c - bv||/||bv||: ", &
!   sqrt(dot_product(r,r)/dot_product(bv,bv))
end function

end module
