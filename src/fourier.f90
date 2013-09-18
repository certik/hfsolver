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

      SUBROUTINE PASSF5_f77(IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
      integer, intent(in) :: IDO, L1
      real(dp), intent(in) :: CC(IDO,5,L1), WA1(:), WA2(:), WA3(:), WA4(:)
      real(dp), intent(out) :: CH(IDO,L1,5)
      real(dp) :: CR2, CR3, CR4, CR5
      real(dp) :: CI2, CI3, CI4, CI5
      real(dp) :: TR2, TR3, TR4, TR5
      real(dp) :: TI2, TI3, TI4, TI5
      real(dp) :: TR11, TI11, TR12, TI12
      real(dp) :: DR2, DR3, DR4, DR5
      real(dp) :: DI2, DI3, DI4, DI5
      integer :: I, K
!     *** TR11=COS(2*PI/5), TI11=-SIN(2*PI/5)
!     *** TR12=-COS(4*PI/5), TI12=-SIN(4*PI/5)
      DATA TR11,TI11,TR12,TI12 /0.3090169943749474241D0, &
           -0.95105651629515357212D0, &
           -0.8090169943749474241D0, -0.58778525229247312917D0/
      IF (IDO .NE. 2) GO TO 102
      DO 101 K=1,L1
         TI5 = CC(2,2,K)-CC(2,5,K)
         TI2 = CC(2,2,K)+CC(2,5,K)
         TI4 = CC(2,3,K)-CC(2,4,K)
         TI3 = CC(2,3,K)+CC(2,4,K)
         TR5 = CC(1,2,K)-CC(1,5,K)
         TR2 = CC(1,2,K)+CC(1,5,K)
         TR4 = CC(1,3,K)-CC(1,4,K)
         TR3 = CC(1,3,K)+CC(1,4,K)
         CH(1,K,1) = CC(1,1,K)+TR2+TR3
         CH(2,K,1) = CC(2,1,K)+TI2+TI3
         CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
         CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
         CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
         CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
         CR5 = TI11*TR5+TI12*TR4
         CI5 = TI11*TI5+TI12*TI4
         CR4 = TI12*TR5-TI11*TR4
         CI4 = TI12*TI5-TI11*TI4
         CH(1,K,2) = CR2-CI5
         CH(1,K,5) = CR2+CI5
         CH(2,K,2) = CI2+CR5
         CH(2,K,3) = CI3+CR4
         CH(1,K,3) = CR3-CI4
         CH(1,K,4) = CR3+CI4
         CH(2,K,4) = CI3-CR4
         CH(2,K,5) = CI2-CR5
  101 CONTINUE
      RETURN
  102 DO 104 K=1,L1
         DO 103 I=2,IDO,2
            TI5 = CC(I,2,K)-CC(I,5,K)
            TI2 = CC(I,2,K)+CC(I,5,K)
            TI4 = CC(I,3,K)-CC(I,4,K)
            TI3 = CC(I,3,K)+CC(I,4,K)
            TR5 = CC(I-1,2,K)-CC(I-1,5,K)
            TR2 = CC(I-1,2,K)+CC(I-1,5,K)
            TR4 = CC(I-1,3,K)-CC(I-1,4,K)
            TR3 = CC(I-1,3,K)+CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
            CH(I,K,1) = CC(I,1,K)+TI2+TI3
            CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
            CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
            CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
            CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
            CR5 = TI11*TR5+TI12*TR4
            CI5 = TI11*TI5+TI12*TI4
            CR4 = TI12*TR5-TI11*TR4
            CI4 = TI12*TI5-TI11*TI4
            DR3 = CR3-CI4
            DR4 = CR3+CI4
            DI3 = CI3+CR4
            DI4 = CI3-CR4
            DR5 = CR2+CI5
            DR2 = CR2-CI5
            DI5 = CI2-CR5
            DI2 = CI2+CR5
            CH(I-1,K,2) = WA1(I-1)*DR2+WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2-WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3+WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3-WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4+WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4-WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5+WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5-WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
      RETURN
      END

      SUBROUTINE PASSF_f77(NAC,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
      integer, intent(out) :: NAC
      integer, intent(in) :: IDO, IP, L1, IDL1
      real(dp), intent(in) :: CC(IDO,IP,L1), WA(:)
      real(dp), intent(out) :: CH(IDO,L1,IP), C1(IDO,L1,IP), C2(IDL1,IP), &
          CH2(IDL1,IP)
      real(dp) :: wai, war
      integer :: I, K, idij, idj, idl, idlj, idot, idp, ik, inc, ipp2, &
          ipph, j, jc, l, lc, nt
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IDP = IP*IDO

      IF (IDO .LT. L1) GO TO 106
      DO 103 J=2,IPPH
         JC = IPP2-J
         DO 102 K=1,L1
            DO 101 I=1,IDO
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  101       CONTINUE
  102    CONTINUE
  103 CONTINUE
      DO 105 K=1,L1
         DO 104 I=1,IDO
            CH(I,K,1) = CC(I,1,K)
  104    CONTINUE
  105 CONTINUE
      GO TO 112
  106 DO 109 J=2,IPPH
         JC = IPP2-J
         DO 108 I=1,IDO
            DO 107 K=1,L1
               CH(I,K,J) = CC(I,J,K)+CC(I,JC,K)
               CH(I,K,JC) = CC(I,J,K)-CC(I,JC,K)
  107       CONTINUE
  108    CONTINUE
  109 CONTINUE
      DO 111 I=1,IDO
         DO 110 K=1,L1
            CH(I,K,1) = CC(I,1,K)
  110    CONTINUE
  111 CONTINUE
  112 IDL = 2-IDO
      INC = 0
      DO 116 L=2,IPPH
         LC = IPP2-L
         IDL = IDL+IDO
         DO 113 IK=1,IDL1
            C2(IK,L) = CH2(IK,1)+WA(IDL-1)*CH2(IK,2)
            C2(IK,LC) = -WA(IDL)*CH2(IK,IP)
  113    CONTINUE
         IDLJ = IDL
         INC = INC+IDO
         DO 115 J=3,IPPH
            JC = IPP2-J
            IDLJ = IDLJ+INC
            IF (IDLJ .GT. IDP) IDLJ = IDLJ-IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            DO 114 IK=1,IDL1
               C2(IK,L) = C2(IK,L)+WAR*CH2(IK,J)
               C2(IK,LC) = C2(IK,LC)-WAI*CH2(IK,JC)
  114       CONTINUE
  115    CONTINUE
  116 CONTINUE
      DO 118 J=2,IPPH
         DO 117 IK=1,IDL1
            CH2(IK,1) = CH2(IK,1)+CH2(IK,J)
  117    CONTINUE
  118 CONTINUE
      DO 120 J=2,IPPH
         JC = IPP2-J
         DO 119 IK=2,IDL1,2
            CH2(IK-1,J) = C2(IK-1,J)-C2(IK,JC)
            CH2(IK-1,JC) = C2(IK-1,J)+C2(IK,JC)
            CH2(IK,J) = C2(IK,J)+C2(IK-1,JC)
            CH2(IK,JC) = C2(IK,J)-C2(IK-1,JC)
  119    CONTINUE
  120 CONTINUE
      NAC = 1
      IF (IDO .EQ. 2) RETURN
      NAC = 0
      DO 121 IK=1,IDL1
         C2(IK,1) = CH2(IK,1)
  121 CONTINUE
      DO 123 J=2,IP
         DO 122 K=1,L1
            C1(1,K,J) = CH(1,K,J)
            C1(2,K,J) = CH(2,K,J)
  122    CONTINUE
  123 CONTINUE
      IF (IDOT .GT. L1) GO TO 127
      IDIJ = 0
      DO 126 J=2,IP
         IDIJ = IDIJ+2
         DO 125 I=4,IDO,2
            IDIJ = IDIJ+2
            DO 124 K=1,L1
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  124       CONTINUE
  125    CONTINUE
  126 CONTINUE
      RETURN
  127 IDJ = 2-IDO
      DO 130 J=2,IP
         IDJ = IDJ+IDO
         DO 129 K=1,L1
            IDIJ = IDJ
            DO 128 I=4,IDO,2
               IDIJ = IDIJ+2
               C1(I-1,K,J) = WA(IDIJ-1)*CH(I-1,K,J)+WA(IDIJ)*CH(I,K,J)
               C1(I,K,J) = WA(IDIJ-1)*CH(I,K,J)-WA(IDIJ)*CH(I-1,K,J)
  128       CONTINUE
  129    CONTINUE
  130 CONTINUE
      RETURN
      END

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
        CH(I,K,1) = &
            (CC(I,1,K) + CC(I,3,K) +    (CC(I,2,K) + CC(I,4,K)))
        CH(I,K,2) = conjg(WA1(I)) * &
            (CC(I,1,K) - CC(I,3,K) - i_*(CC(I,2,K) - CC(I,4,K)))
        CH(I,K,3) = conjg(WA2(I)) * &
            (CC(I,1,K) + CC(I,3,K) -    (CC(I,2,K) + CC(I,4,K)))
        CH(I,K,4) = conjg(WA3(I)) * &
            (CC(I,1,K) - CC(I,3,K) + i_*(CC(I,2,K) - CC(I,4,K)))
    end do
end do
end subroutine

subroutine passf5(IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
! FFT pass of factor 2
use iso_c_binding, only: c_f_pointer, c_loc
integer, intent(in) :: IDO, L1
complex(dp), intent(in), target :: CC(IDO, 5, L1), WA1(:), WA2(:), WA3(:), &
    WA4(:)
complex(dp), intent(out), target :: CH(IDO, L1, 5)
complex(dp), target :: WA1p(size(WA1)), WA2p(size(WA2)), WA3p(size(WA3)), &
    WA4p(size(WA4))
real(dp), pointer :: CC_r(:, :, :), CH_r(:, :, :), WA1_r(:), WA2_r(:), &
    WA3_r(:), WA4_r(:)
call c_f_pointer(c_loc(CC), CC_r, [size(CC)*2])
call c_f_pointer(c_loc(CH), CH_r, [size(CH)*2])
WA1p = WA1
WA2p = WA2
WA3p = WA3
WA4p = WA4
call c_f_pointer(c_loc(WA1p), WA1_r, [size(WA1)*2])
call c_f_pointer(c_loc(WA2p), WA2_r, [size(WA2)*2])
call c_f_pointer(c_loc(WA3p), WA3_r, [size(WA3)*2])
call c_f_pointer(c_loc(WA4p), WA4_r, [size(WA4)*2])
call PASSF5_f77(IDO*2, L1, CC_r, CH_r, WA1_r, WA2_r, WA3_r, WA4_r)
end subroutine

subroutine passf(NAC, IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
use iso_c_binding, only: c_f_pointer, c_loc
integer, intent(out) :: NAC
integer, intent(in) :: IDO, IP, L1, IDL1
complex(dp), intent(inout), target :: CC(IDO,IP,L1)
complex(dp), intent(in), target ::  WA(:)
complex(dp), intent(out), target :: CH(IDO,L1,IP), C1(IDO,L1,IP), C2(IDL1,IP), &
          CH2(IDL1,IP)
complex(dp), target :: WA1p(size(WA))
real(dp), pointer :: CC_r(:, :, :), CH_r(:, :, :), C1_r(:, :, :), C2_r(:, :), &
    CH2_r(:, :), WA1_r(:)
call c_f_pointer(c_loc(CC), CC_r, [size(CC)*2])
call c_f_pointer(c_loc(CH), CH_r, [size(CH)*2])
call c_f_pointer(c_loc(C1), C1_r, [size(C1)*2])
call c_f_pointer(c_loc(C2), C2_r, [size(C2)*2])
call c_f_pointer(c_loc(CH2), CH2_r, [size(CH2)*2])
WA1p = WA
call c_f_pointer(c_loc(WA1p), WA1_r, [size(WA)*2])
call PASSF_f77(NAC, IDO*2, IP, L1, IDL1, CC_r, C1_r, C2_r, CH_r, CH2_r, WA1_r)
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


SUBROUTINE CFFTF1(N,C,CH,WA,IFAC)
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
            call passf2(IDO,L1,C,CH,conjg(w(:, 1)))
        else
            call passf2(IDO,L1,CH,C,conjg(w(:, 1)))
        end if
    case (3)
        if (NA == 0) then
            call passf3(IDO,L1,C,CH,conjg(w(:, 1)),conjg(w(:, 2)))
        else
            call passf3(IDO,L1,CH,C,conjg(w(:, 1)),conjg(w(:, 2)))
        end if
    case (5)
        if (NA == 0) then
            call passf5(IDO,L1,C,CH,w(:, 1),w(:, 2),w(:, 3),w(:, 4))
        else
            call passf5(IDO,L1,CH,C,w(:, 1),w(:, 2),w(:, 3),w(:, 4))
        end if
    case default
        IDL1 = 2*IDO*L1
        if (NA == 0) then
            call passf(NAC, IDO, IP, L1, IDL1, C, C, C, CH, CH, &
                WA(IW:IW+(IP-1)*IDO-1))
        else
            call passf(NAC, IDO, IP, L1, IDL1, CH, CH, CH, C, C, &
                WA(IW:IW+(IP-1)*IDO-1))
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
