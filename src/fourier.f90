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
integer :: n, i, j, k1, ii, i1, ido, idot, ipm, ip, l1, l2, ld
real(dp) :: argh, arg, argld, fi
n = size(wa)

      ARGH = 2*pi/n
      I = 2
      L1 = 1
      do K1 = 1, size(fac)
         IP = fac(K1)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO+IDO+2
         IPM = IP-1
         DO 109 J=1,IPM
            I1 = I
            wa(i/2) = 1
            LD = LD+L1
            FI = 0
            ARGLD = LD*ARGH
            DO II=4,IDOT,2
               I = I+2
               FI = FI+1
               ARG = FI*ARGLD
               wa(i/2) = cos(arg) + i_ * sin(arg)
            end do
            IF (IP .LE. 5) GO TO 109
            wa(i1/2) = wa(i/2)
  109    CONTINUE
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

      SUBROUTINE PASSF4_f77(IDO,L1,CC,CH,WA1,WA2,WA3)
      integer, intent(in) :: IDO, L1
      real(dp), intent(in) :: CC(IDO,4,L1), WA1(:), WA2(:), WA3(:)
      real(dp), intent(out) :: CH(IDO,L1,4)
      real(dp) :: CR2, CR3, CR4
      real(dp) :: CI2, CI3, CI4
      real(dp) :: TR1, TR2, TR3, TR4
      real(dp) :: TI1, TI2, TI3, TI4
      integer :: I, K
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI1 = CC(2,1,K)-CC(2,3,K)
         TI2 = CC(2,1,K)+CC(2,3,K)
         TR4 = CC(2,2,K)-CC(2,4,K)
         TI3 = CC(2,2,K)+CC(2,4,K)
         TR1 = CC(1,1,K)-CC(1,3,K)
         TR2 = CC(1,1,K)+CC(1,3,K)
         TI4 = CC(1,4,K)-CC(1,2,K)
         TR3 = CC(1,2,K)+CC(1,4,K)
         CH(1,K,1) = TR2+TR3
         CH(1,K,3) = TR2-TR3
         CH(2,K,1) = TI2+TI3
         CH(2,K,3) = TI2-TI3
         CH(1,K,2) = TR1+TR4
         CH(1,K,4) = TR1-TR4
         CH(2,K,2) = TI1+TI4
         CH(2,K,4) = TI1-TI4
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI1 = CC(I,1,K)-CC(I,3,K)
            TI2 = CC(I,1,K)+CC(I,3,K)
            TI3 = CC(I,2,K)+CC(I,4,K)
            TR4 = CC(I,2,K)-CC(I,4,K)
            TR1 = CC(I-1,1,K)-CC(I-1,3,K)
            TR2 = CC(I-1,1,K)+CC(I-1,3,K)
            TI4 = CC(I-1,4,K)-CC(I-1,2,K)
            TR3 = CC(I-1,2,K)+CC(I-1,4,K)
            CH(I-1,K,1) = TR2+TR3
            CR3 = TR2-TR3
            CH(I,K,1) = TI2+TI3
            CI3 = TI2-TI3
            CR2 = TR1+TR4
            CR4 = TR1-TR4
            CI2 = TI1+TI4
            CI4 = TI1-TI4
            CH(I-1,K,2) = WA1(I-1)*CR2+WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2-WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3+WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3-WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4+WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4-WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

subroutine passf4(IDO, L1, CC, CH, WA1, WA2, WA3)
! FFT pass of factor 2
use iso_c_binding, only: c_f_pointer, c_loc
integer, intent(in) :: IDO, L1
complex(dp), intent(in), target :: CC(IDO, 4, L1), WA1(:), WA2(:), WA3(:)
complex(dp), intent(out), target :: CH(IDO, L1, 4)
complex(dp), target :: WA1p(size(WA1)), WA2p(size(WA2)), WA3p(size(WA3))
real(dp), pointer :: CC_r(:, :, :), CH_r(:, :, :), WA1_r(:), WA2_r(:), WA3_r(:)
call c_f_pointer(c_loc(CC), CC_r, [size(CC)*2])
call c_f_pointer(c_loc(CH), CH_r, [size(CH)*2])
WA1p = WA1
WA2p = WA2
WA3p = WA3
call c_f_pointer(c_loc(WA1p), WA1_r, [size(WA1)*2])
call c_f_pointer(c_loc(WA2p), WA2_r, [size(WA2)*2])
call c_f_pointer(c_loc(WA3p), WA3_r, [size(WA3)*2])
call PASSF4_f77(IDO*2, L1, CC_r, CH_r, WA1_r, WA2_r, WA3_r)
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


      SUBROUTINE CFFTF1 (N,C,CH,WA,IFAC)
      ! NOTE: PASSF3, PASSF5 and PASSF are commented out for now
      integer, intent(in) :: N
      complex(dp), intent(inout) :: C(:)
      complex(dp), intent(out) :: CH(:)
      complex(dp), intent(in) :: WA(:)
      integer, intent(in) :: IFAC(:)
      integer :: k1, l1, na, iw, ip, l2, ido, idot, idl1, ix2, ix3, ix4, &
          n2
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,size(ifac)
         IP = IFAC(K1)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO+IDO
         IDL1 = IDOT*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IF (NA .NE. 0) GO TO 101
         call passf4(IDOT/2,L1,C,CH,WA(IW/2+1:),WA(IX2/2+1:),WA(IX3/2+1:))
         GO TO 102
  101    call passf4(IDOT/2,L1,CH,C,WA(IW/2+1:),WA(IX2/2+1:),WA(IX3/2+1:))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         call passf2(IDOT/2,L1,C,CH,conjg(WA(IW/2+1:)))
         GO TO 105
  104    call passf2(IDOT/2,L1,CH,C,conjg(WA(IW/2+1:)))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDOT
         IF (NA .NE. 0) GO TO 107
         !CALL PASSF3 (IDOT,L1,C,CH,WA(IW),WA(IX2))
         GO TO 108
  107    continue !CALL PASSF3 (IDOT,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDOT
         IX3 = IX2+IDOT
         IX4 = IX3+IDOT
         IF (NA .NE. 0) GO TO 110
         !CALL PASSF5 (IDOT,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110    continue !CALL PASSF5 (IDOT,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         !CALL PASSF (NAC,IDOT,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
         GO TO 114
  113    continue !CALL PASSF (NAC,IDOT,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    continue !IF (NAC .NE. 0) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDOT
  116 CONTINUE
      IF (NA .EQ. 0) RETURN
      N2 = N+N
      C = CH
      END

function fft_pass(x) result(p)
real(dp), intent(in) :: x(:)
complex(dp) :: p(size(x))
p = x
call fft_pass_inplace(p)
end function

end module
