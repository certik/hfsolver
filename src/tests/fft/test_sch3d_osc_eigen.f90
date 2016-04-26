program test_sch3d_osc_eigen

! Tests 1D harmonic oscillator ground state minimization (imaginary time
! propagation).

use types, only: dp
use constants, only: i_
use ofdft, only: read_pseudo
use ofdft_fft, only: free_energy, radial_potential_fourier, &
    reciprocal_space_vectors, free_energy_min, real2fourier, integral, &
    fourier2real, real_space_vectors, vtk_save, integralG
use utils, only: loadtxt, stop_error, assert, linspace, strfmt
use splines, only: spline3pars, iixmin, poly3, spline3ders
use interp3d, only: trilinear
use md, only: positions_fcc, positions_bcc
use linalg, only: eigvals
use sorting, only: argsort
implicit none
integer :: Ng
real(dp), allocatable :: G(:,:,:,:), G2(:,:,:)
real(dp), allocatable :: Xn(:,:,:,:), Vn(:,:,:), r(:,:,:)
complex(dp), allocatable, dimension(:) :: lam
complex(dp), allocatable :: H(:,:), psi(:,:,:)
integer, allocatable :: idx(:)
real(dp) :: L
real(dp) :: omega
integer :: i, ai, aj, ak, bi, bj, bk, ci, cj, ck, a_idx, b_idx !, a, b

Ng = 8

L = 10._dp

allocate(G(Ng,Ng,Ng,3), G2(Ng,Ng,Ng), psi(Ng,Ng,Ng))
allocate(Xn(Ng, Ng, Ng, 3), Vn(Ng, Ng, Ng), r(Ng, Ng, Ng))
allocate(H(Ng**3, Ng**3))
allocate(lam(Ng**3), idx(Ng**3))

call real_space_vectors(L, Xn)
call reciprocal_space_vectors(L, G, G2)
omega = 1._dp
r = sqrt(sum((Xn-L/2)**2, dim=4))
Vn = omega**2 * r**2 / 2

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
print *, "Eigenvalues:"
do i = Ng**3, 1, -1
    print *, i, lam(idx(i))
end do

print *, "E_tot_exact =", 3*omega/2

end program
