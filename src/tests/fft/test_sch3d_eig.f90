program test_sch3d_eig

use types, only: dp
use constants, only: i_
use ofdft, only: read_pseudo
use ofdft_fft, only: free_energy, radial_potential_fourier, &
    reciprocal_space_vectors, free_energy_min, real2fourier, integral, &
    fourier2real, real_space_vectors, vtk_save, integralG
use utils, only: loadtxt, stop_error, assert, linspace, strfmt, init_random
use splines, only: spline3pars, iixmin, poly3, spline3ders
use interp3d, only: trilinear
use md, only: positions_fcc, positions_bcc
use linalg, only: eigvals, inv, solve
use arpack, only: eig
use sorting, only: argsort
implicit none
integer :: Ng
real(dp), allocatable :: G(:,:,:,:), G2(:,:,:)
real(dp), allocatable :: Xn(:,:,:,:), Vn(:,:,:), r(:,:,:)
complex(dp), allocatable, dimension(:) :: lam
complex(dp), allocatable :: H(:,:), psi(:,:,:), psiG(:,:,:)
integer, allocatable :: idx(:)
real(dp) :: L
real(dp) :: lambda !, omega
integer :: i, ai, aj, ak, bi, bj, bk, ci, cj, ck, a_idx, b_idx !, a, b
complex(dp) :: lam0
integer :: n, nev, ncv
real(dp), allocatable :: d(:), v(:,:)

Ng = 8

L = 10._dp

allocate(G(Ng,Ng,Ng,3), G2(Ng,Ng,Ng), psi(Ng,Ng,Ng), psiG(Ng,Ng,Ng))
allocate(Xn(Ng, Ng, Ng, 3), Vn(Ng, Ng, Ng), r(Ng, Ng, Ng))
allocate(H(Ng**3, Ng**3))
allocate(lam(Ng**3), idx(Ng**3))

call real_space_vectors([L, L, L], Xn)
call reciprocal_space_vectors([L, L, L], G, G2)
!omega = 1._dp
r = sqrt(sum((Xn-L/2)**2, dim=4))
!Vn = omega**2 * r**2 / 2
lambda = 0.2_dp
Vn = -exp(-lambda*r**2)

!a = 2
!b = 12
!Vn = -1/sqrt(a+(Xn-L/2)**b)


call fourier2real(G2/2+0*i_, psi)

!do j = 1, Ng
!do i = 1, Ng
!    k = i-j+1
!    if (k < 1) k = k + Ng
!    H(i,j) = psi(k) / Ng
!end do
!end do

print *, "Assembly"

a_idx = 0
do ak = 1, Ng
do aj = 1, Ng
do ai = 1, Ng
    a_idx = a_idx+1
    b_idx = 0
    do bk = 1, Ng
    do bj = 1, Ng
    do bi = 1, Ng
        b_idx = b_idx+1
        ci = ai-bi+1
        cj = aj-bj+1
        ck = ak-bk+1
        if (ci < 1) ci = ci + Ng
        if (cj < 1) cj = cj + Ng
        if (ck < 1) ck = ck + Ng
        H(a_idx,b_idx) = psi(ci,cj,ck) / Ng**3
    end do
    end do
    end do
end do
end do
end do

a_idx = 0
do ak = 1, Ng
do aj = 1, Ng
do ai = 1, Ng
    a_idx = a_idx+1
    H(a_idx,a_idx) = H(a_idx,a_idx) + Vn(ai,aj,ak)
end do
end do
end do

print *, "Eigensolver"
lam = eigvals(H)

print *, "Sorting"
idx = argsort(real(lam, dp))

print *, "Power iteration"
call init_random()
lam0 = power_iteration(H)
print *, lam0
print *, "Eigenvalues:"
do i = Ng**3, 1, -1
    print *, i, lam(idx(i))
end do

print *, "Inverse iteration"
lam0 = inverse_iteration(H, -0.20_dp)
print *, lam0

print *, "Arpack"

n = Ng**3
nev = 4
ncv = 100
allocate(v(n,ncv), d(ncv))
call eig(n, nev, ncv, "SA", av, d, v)
do i = 1, nev
    print *, d(i)
end do

!print *, "E_tot_exact =", 3*omega/2

contains

    complex(dp) function power_iteration(A) result(lam)
    complex(dp), intent(in) :: A(:, :)
    complex(dp) :: y(size(A, 1)), v(size(A, 1))
    real(dp) :: y0(size(A, 1))
    integer :: i
    call random_number(y0)
    y = y0
    print *
    do i = 1, 40
        v = y / sqrt(sum(abs(y)**2))
        y = matmul(A, v)
        lam = dot_product(conjg(v), y)
        print *, i, lam
    end do
    end function

    complex(dp) function inverse_iteration(A, mu) result(lam)
    complex(dp), intent(in) :: A(:, :)
    real(dp), intent(in) :: mu
    complex(dp) :: M(size(A, 1), size(A, 2))
    complex(dp) :: y(size(A, 1)), v(size(A, 1)), theta
    real(dp) :: y0(size(A, 1))
    integer :: i, N
    logical, parameter :: invert = .false.
    N = size(A, 1)
    call assert(size(A, 2) == N)
    M = A
    do i = 1, N
        M(i,i) = M(i,i) - mu
    end do
    if (invert) then
        print *, "inverting:"
        M = inv(M)
        print *, "done"
    end if

    call random_number(y0)
    y = y0
    print *
    do i = 1, 15
        v = y / sqrt(sum(abs(y)**2))
        if (invert) then
            y = matmul(M, v)
        else
            y = solve(M, v)
        end if
        theta = dot_product(conjg(v), y)
        lam = mu + 1/theta
        print *, i, lam
    end do
    end function

    subroutine av(x, y)
    ! Compute y = A*x
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    call real2fourier(reshape(x, [Ng,Ng,Ng]), psiG)
    call fourier2real(G2/2*psiG, psi)
    !print *, maxval(abs(aimag(psi)))
    y = reshape(real(psi,dp), [Ng**3]) + reshape(Vn, [Ng**3])*x
    end

end program
