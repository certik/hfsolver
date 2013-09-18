module fourier

! Calculates the discrete Fourier transform using the simple direct but slow
! O(n^2) algorithm.

use types, only: dp
use constants, only: i_, pi
use utils, only: stop_error, assert
implicit none
private
public dft, idft, fft, fft_vectorized, fft_pass, fft_pass_inplace, &
    fft_vectorized_inplace, calculate_factors

contains

function dft(x) result(p)
! Compute the one-dimensional discrete Fourier transform
real(dp), intent(in) :: x(:)
complex(dp) :: p(size(x))
complex(dp) :: F(size(x), size(x))
integer :: N, i, j
N = size(x)
forall(i=0:N-1, j=0:N-1, i >= j) F(i+1, j+1) = exp(-2*pi*i_*i*j / N)
forall(i=1:N, j=1:N, i < j) F(i, j) = F(j, i)
p = matmul(F, x)
end function

function idft(p) result(x)
! Compute the one-dimensional inverse discrete Fourier transform
! The normalization is such that idft(dft(x)) == x to within numerical
! accuracy.
complex(dp), intent(in) :: p(:)
complex(dp) :: x(size(p))
complex(dp) :: F(size(p), size(p))
integer :: N, i, j
N = size(p)
forall(i=0:N-1, j=0:N-1, i >= j) F(i+1, j+1) = exp(+2*pi*i_*i*j/ N)
forall(i=1:N, j=1:N, i < j) F(i, j) = F(j, i)
x = matmul(F, p) / N
end function

recursive function fft(x) result(p)
! A recursive implementation of the 1D Cooley-Tukey FFT
real(dp), intent(in) :: x(:)
complex(dp) :: p(size(x)), X_even(size(x)/2), X_odd(size(x)/2)
complex(dp) :: factor(size(x))
integer :: N, i
N = size(x)
if (iand(N, N-1) /= 0) call stop_error("size of x must be a power of 2")
if (N <= 4) then
    p = dft(x)
else
    X_even = fft(x(::2))
    X_odd = fft(x(2::2))
    forall(i=0:N-1) factor(i+1) = exp(-2*pi*i_*i/N)
    p = [X_even + factor(:N / 2) * X_odd, X_even + factor(N / 2+1:) * X_odd]
end if
end function

subroutine dft_vec(Ns, N, x, p)
! Compute the one-dimensional discrete Fourier transform on each row of 'x'
! separately.
integer, intent(in) :: Ns, N
complex(dp), intent(in) :: x(Ns, N)
complex(dp), intent(out) :: p(Ns, N)
complex(dp) :: F(N, N)
integer :: i, j
forall(i=0:N-1, j=0:N-1, i >= j) F(i+1, j+1) = exp(-2*pi*i_*i*j / N)
forall(i=1:N, j=1:N, i < j) F(i, j) = F(j, i)
p = matmul(x, F)
end subroutine

subroutine fft_step(Ns, Nmin, x, p)
integer, intent(in) :: Ns, Nmin
complex(dp), intent(in) :: x(Ns, Nmin)
complex(dp), intent(out) :: p(Ns/2, Nmin*2)
complex(dp) :: tmp
integer :: i
do i = 1, Nmin
    !tmp = exp(-pi*i_*(i-1)/Nmin)
    ! The same as the previous line, just faster:
    tmp = cos(pi*(i-1)/Nmin) - i_*sin(pi*(i-1)/Nmin)
    p(:,      i) = x(:Ns/2, i) + tmp * x(Ns/2+1:, i)
    p(:, Nmin+i) = x(:Ns/2, i) - tmp * x(Ns/2+1:, i)
end do
end subroutine

subroutine fft_vectorized_inplace(x)
! A vectorized, non-recursive version of the Cooley-Tukey FFT
complex(dp), intent(inout), target :: x(:)
complex(dp), target :: p(size(x))
integer :: N, Nmin, Ns
logical :: p_is_result
N = size(x)
if (iand(N, N-1) /= 0) call stop_error("size of x must be a power of 2")
Nmin = 1
Ns = N
p_is_result = .true.
do while (Nmin < N)
    if (p_is_result) then
        call fft_step(Ns, Nmin, x, p)
    else
        call fft_step(Ns, Nmin, p, x)
    end if
    Nmin = Nmin * 2
    Ns = Ns / 2
    p_is_result = .not. p_is_result
end do
if (.not. p_is_result) x = p
end subroutine

function fft_vectorized(x) result(p)
real(dp), intent(in) :: x(:)
complex(dp) :: p(size(x))
p = x
call fft_vectorized_inplace(p)
end function

subroutine precalculate_coeffs(wa)
! Precalculates all cos/sin factors
complex(dp), intent(out) :: wa(:)
integer :: n, k, i, idx
n = size(wa)
k = n / 2
idx = 1
do while (k > 0)
    wa(idx) = 1
    do i = 1, k
        idx = idx + 1
        ! Equivalent to exp(-i_*i*pi/k) but faster:
        wa(idx) = cos(i*pi/k) - i_ * sin(i*pi/k)
    end do
    k = k/2
end do
end subroutine

subroutine precalculate_angles(fac, wa)
! Precalculates all cos/sin factors
integer, intent(in) :: fac(:)
complex(dp), intent(out) :: wa(:)
integer :: fi, n, i, j, k1, i1, ip, l1, l2, ld
real(dp) :: arg, argld
n = size(wa)
I = 1
L1 = 1
do K1 = 1, size(fac)
    IP = fac(K1)
    LD = 0
    L2 = L1*IP
    do J = 1, IP-1
        I1 = I
        wa(i) = 1
        LD = LD+L1
        ARGLD = 2*pi*LD/n
        do FI = 1, N/L2
           I = I+1
           ARG = FI*ARGLD
           wa(i) = cos(arg) - i_ * sin(arg)
        end do
        if (IP > 5) wa(i1) = wa(i)
    end do
    L1 = L2
end do
end subroutine

subroutine passf2(IDO, L1, CC, CH, WA1)
! FFT pass of factor 2
integer, intent(in) :: IDO, L1
complex(dp), intent(in) :: CC(IDO, 2, L1), WA1(:)
complex(dp), intent(out) :: CH(IDO, L1, 2)
integer :: I, K
do K = 1, L1
    do I = 1, IDO
        CH(I,K,1) = CC(I,1,K)+CC(I,2,K)
        CH(I,K,2) = WA1(I)*(CC(I,1,K)-CC(I,2,K))
    end do
end do
end subroutine

subroutine passf3(IDO, L1, CC, CH, WA1, WA2)
integer, intent(in) :: IDO, L1
complex(dp), intent(in) :: CC(IDO,3,L1), WA1(:), WA2(:)
complex(dp), intent(out) :: CH(IDO,L1,3)
complex(dp) :: C1, C2, C3
integer :: I, K
real(dp), parameter :: taur = -0.5_dp, taui = -sqrt(3._dp)/2
do K=1,L1
    do I=1,IDO
        CH(I,K,1) = CC(I,1,K)+CC(I,2,K)+CC(I,3,K)
        C1 = CC(I,1,K) + TAUR * (CC(I,2,K)+CC(I,3,K))
        C3 = TAUI * (CC(I,2,K)-CC(I,3,K))
        ! The same as C2 = C3 * exp(i_*pi/2) but faster:
        C2 = - aimag(C3) + i_ * real(C3, dp)
        CH(I,K,2) = WA1(I) * (C1 + C2)
        CH(I,K,3) = WA2(I) * (C1 - C2)
    end do
end do
end subroutine

subroutine passf4(IDO, L1, CC, CH, WA1, WA2, WA3)
integer, intent(in) :: IDO, L1
complex(dp), intent(in) :: CC(IDO, 4, L1), WA1(:), WA2(:), WA3(:)
complex(dp), intent(out) :: CH(IDO, L1, 4)
integer :: I, K
do K = 1, L1
    do I = 1, IDO
        CH(I,K,1) =          (CC(I,1,K) + CC(I,3,K) +   (CC(I,2,K) + CC(I,4,K)))
        CH(I,K,2) = WA1(I) * (CC(I,1,K) - CC(I,3,K) -i_*(CC(I,2,K) - CC(I,4,K)))
        CH(I,K,3) = WA2(I) * (CC(I,1,K) + CC(I,3,K) -   (CC(I,2,K) + CC(I,4,K)))
        CH(I,K,4) = WA3(I) * (CC(I,1,K) - CC(I,3,K) +i_*(CC(I,2,K) - CC(I,4,K)))
    end do
end do
end subroutine

subroutine passf5(IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
integer, intent(in) :: IDO, L1
complex(dp), intent(in) :: CC(IDO, 5, L1), WA1(:), WA2(:), WA3(:), WA4(:)
complex(dp), intent(out) :: CH(IDO, L1, 5)
complex(dp) :: C2, C3, C4, C5, C4b, C5b
integer :: I, K
real(dp), parameter :: tr11 =  cos(2*pi/5)
real(dp), parameter :: ti11 = -sin(2*pi/5)
real(dp), parameter :: tr12 =  cos(4*pi/5)
real(dp), parameter :: ti12 = -sin(4*pi/5)
do K = 1, L1
    do I = 1, IDO
        CH(I,K,1) = CC(I,1,K) + CC(I,2,K) + CC(I,3,K) + CC(I,4,K) + CC(I,5,K)
        C2 = CC(I,1,K) + TR11*(CC(I,2,K)+CC(I,5,K)) + TR12*(CC(I,3,K)+CC(I,4,K))
        C3 = CC(I,1,K) + TR12*(CC(I,2,K)+CC(I,5,K)) + TR11*(CC(I,3,K)+CC(I,4,K))
        C4 = TI12*(CC(I,2,K)-CC(I,5,K)) - TI11*(CC(I,3,K)-CC(I,4,K))
        C5 = TI11*(CC(I,2,K)-CC(I,5,K)) + TI12*(CC(I,3,K)-CC(I,4,K))
        C4b = aimag(C4) - i_*real(C4, dp)
        C5b = aimag(C5) - i_*real(C5, dp)
        CH(I,K,2) = WA1(I) * (C2 - C5b)
        CH(I,K,3) = WA2(I) * (C3 - C4b)
        CH(I,K,4) = WA3(I) * (C3 + C4b)
        CH(I,K,5) = WA4(I) * (C2 + C5b)
    end do
end do
end subroutine

subroutine passf(NAC, IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
! Works for any odd IP
integer, intent(out) :: NAC
integer, intent(in) :: IDO, IP, L1, IDL1
complex(dp), intent(inout) :: CC(IDO,IP,L1)
complex(dp), intent(in) ::  WA(:)
complex(dp), intent(out) :: CH(IDO,L1,IP), C1(IDO,L1,IP), C2(IDL1,IP), &
          CH2(IDL1,IP)
integer :: I, K, idij, idj, idl, idlj, idp, ik, inc, ipp2, &
    ipph, j, jc, l, lc, nt
NT = IP*IDL1
IPP2 = IP+2
IPPH = (IP+1)/2
IDP = IP*IDO
do J = 2, IPPH
    JC = IPP2-J
    CH(:,:,J ) = CC(:,J,:)+CC(:,JC,:)
    CH(:,:,JC) = CC(:,J,:)-CC(:,JC,:)
end do
CH(:,:,1) = CC(:,1,:)
IDL = 1-IDO
INC = 0
do L=2,IPPH
    LC = IPP2-L
    IDL = IDL+IDO
    do IK=1,IDL1
        C2(IK,L ) = CH2(IK,1) +  real(WA(IDL),dp) * CH2(IK,2 )
        C2(IK,LC) =           - aimag(WA(IDL))    * CH2(IK,IP)
    end do
    IDLJ = IDL
    INC = INC+IDO
    do J=3,IPPH
        JC = IPP2-J
        IDLJ = IDLJ+INC
        if (IDLJ > IDP) IDLJ = IDLJ-IDP
        do IK=1,IDL1
            C2(IK,L ) = C2(IK,L ) +  real(WA(IDLJ),dp) * CH2(IK,J )
            C2(IK,LC) = C2(IK,LC) - aimag(WA(IDLJ))    * CH2(IK,JC)
        end do
    end do
end do
do J = 2, IPPH
    CH2(:,1) = CH2(:,1)+CH2(:,J)
end do
do J=2,IPPH
    JC = IPP2-J
    do IK=1,IDL1
        CH2(IK,J ) = C2(IK,J) - aimag(C2(IK,JC)) + i_ * real(C2(IK,JC))
        CH2(IK,JC) = C2(IK,J) + aimag(C2(IK,JC)) - i_ * real(C2(IK,JC))
    end do
end do
NAC = 1
if (IDO == 1) return
NAC = 0
C2(:,1) = CH2(:,1)
C1(1,:,:) = CH(1,:,:)
IDJ = 1-IDO
do J=2,IP
    IDJ = IDJ+IDO
    do K=1,L1
        IDIJ = IDJ
        do I=2,IDO
            IDIJ = IDIJ+1
            C1(I,K,J) = conjg(WA(IDIJ))*CH(I,K,J)
        end do
    end do
end do
end subroutine

subroutine calculate_factors(n, fac)
integer, intent(in) :: n
integer, intent(out), allocatable :: fac(:)
! TODO: add checks that we don't go over MAX_LENGTH below:
integer, parameter :: MAX_LENGTH = 1000
integer :: fac_tmp(MAX_LENGTH)
integer, parameter :: NTRYH(*) = [3, 4, 2, 5]
integer :: NL, NF, I, J, NTRY, IB, NQ, NR
if (n == 1) then
    allocate(fac(1))
    fac(1) = 1
    return
end if
NL = N
NF = 0
J = 0
do while (NL /= 1)
    J = J+1
    IF (J <= 4) then
        NTRY = NTRYH(J)
    else
        NTRY = NTRY+2
    end if
    ! Divide by NTRY as many times as we can:
    do while (NL /= 1)
        NQ = NL/NTRY
        NR = NL-NTRY*NQ
        IF (NR /= 0) exit
        NF = NF+1
        fac_tmp(NF) = NTRY
        NL = NQ
        if (NTRY == 2 .and. NF > 1) then
            do I = 2, NF
                IB = NF-I+2
                fac_tmp(IB) = fac_tmp(IB-1)
            end do
            fac_tmp(1) = 2
        end if
    end do
end do
allocate(fac(NF))
fac = fac_tmp(:NF)
end subroutine

subroutine fft_pass_inplace(x)
complex(dp), intent(inout) :: x(:)
complex(dp), dimension(size(x)) :: angles, CH
integer, allocatable :: fac(:)
integer :: n
n = size(x)
call calculate_factors(n, fac)
call precalculate_angles(fac, angles)
call cfftf1(n, x, CH, angles, fac)
end subroutine


subroutine cfftf1(N,C,CH,WA,IFAC)
integer, intent(in) :: N
complex(dp), intent(inout) :: C(:)
complex(dp), intent(out) :: CH(:)
complex(dp), intent(in), target :: WA(:)
integer, intent(in) :: IFAC(:)
integer :: k1, l1, na, iw, ip, l2, ido, idl1, nac
complex(dp), pointer :: w(:, :)
NA = 0
L1 = 1
IW = 1
do K1 = 1, size(ifac)
    IP = IFAC(K1)
    L2 = IP*L1
    IDO = N/L2
    w(1:IDO,1:IP-1) => WA(IW:IW+(IP-1)*IDO-1)
    select case(IP)
    case (4)
        if (NA == 0) then
            call passf4(IDO,L1,C,CH,w(:, 1),w(:, 2),w(:, 3))
        else
            call passf4(IDO,L1,CH,C,w(:, 1),w(:, 2),w(:, 3))
        end if
    case (2)
        if (NA == 0) then
            call passf2(IDO,L1,C,CH,w(:, 1))
        else
            call passf2(IDO,L1,CH,C,w(:, 1))
        end if
    case (3)
        if (NA == 0) then
            call passf3(IDO,L1,C,CH,w(:, 1),w(:, 2))
        else
            call passf3(IDO,L1,CH,C,w(:, 1),w(:, 2))
        end if
    case (5)
        if (NA == 0) then
            call passf5(IDO,L1,C,CH,w(:, 1),w(:, 2),w(:, 3),w(:, 4))
        else
            call passf5(IDO,L1,CH,C,w(:, 1),w(:, 2),w(:, 3),w(:, 4))
        end if
    case default
        IDL1 = IDO*L1
        if (NA == 0) then
            call passf(NAC, IDO, IP, L1, IDL1, C, C, C, CH, CH, &
                conjg(WA(IW:IW+(IP-1)*IDO-1)))
        else
            call passf(NAC, IDO, IP, L1, IDL1, CH, CH, CH, C, C, &
                conjg(WA(IW:IW+(IP-1)*IDO-1)))
        end if
        IF (NAC == 0) NA = 1-NA
    end select
    NA = 1-NA
    L1 = L2
    IW = IW+(IP-1)*IDO
end do
if (NA /= 0) C = CH
end subroutine

function fft_pass(x) result(p)
real(dp), intent(in) :: x(:)
complex(dp) :: p(size(x))
p = x
call fft_pass_inplace(p)
end function

end module
