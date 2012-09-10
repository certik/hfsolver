module debye_potential_series
use types, only: dp
implicit none

integer, parameter :: qp = 8

! All these assume that r1 >= r2

contains

!******************************************************************************
!*                    Code generated with sympy 0.7.1-git                     *
!*                                                                            *
!*              See http://www.sympy.org/ for more information.               *
!*                                                                            *
!*                       This file is part of 'project'                       *
!******************************************************************************

real(qp) function S0(alpha, t)
implicit none
real(qp), intent(in) :: alpha
real(qp), intent(in) :: t

S0 = alpha**12*t**12/6227020800._qp + alpha**10*t**10/39916800 + alpha**8*t &
      **8/362880 + alpha**6*t**6/5040 + alpha**4*t**4/120 + alpha**2*t &
      **2/6 + 1

end function

real(qp) function S1(alpha, t)
implicit none
real(qp), intent(in) :: alpha
real(qp), intent(in) :: t

S1 = alpha**13*t**13/31135104000._qp + alpha**12*t**13/31135104000._qp + alpha** &
      11*t**11/172972800 + alpha**10*t**11/172972800 + alpha**9*t**9/ &
      1330560 + alpha**8*t**9/1330560 + alpha**7*t**7/15120 + alpha**6* &
      t**7/15120 + alpha**5*t**5/280 + alpha**4*t**5/280 + alpha**3*t** &
      3/10 + alpha**2*t**3/10 + alpha*t + t

end function

real(qp) function S2(alpha, t)
implicit none
real(qp), intent(in) :: alpha
real(qp), intent(in) :: t

S2 = alpha**13*t**14/105859353600._qp + alpha**12*(t**14/105859353600._qp + t** &
      12/1556755200) + alpha**11*t**12/518918400 + alpha**10*(t**12/ &
      518918400 + t**10/10378368) + alpha**9*t**10/3459456 + alpha**8*( &
      t**10/3459456 + t**8/99792) + alpha**7*t**8/33264 + alpha**6*(t** &
      8/33264 + t**6/1512) + alpha**5*t**6/504 + alpha**4*(t**6/504 + t &
      **4/42) + alpha**3*t**4/14 + alpha**2*(t**4/14 + t**2/3) + alpha* &
      t**2 + t**2

end function

real(qp) function S3(alpha, t)
implicit none
real(qp), intent(in) :: alpha
real(qp), intent(in) :: t

S3 = alpha**13*(t**15/287332531200._qp + t**13/18903456000._qp) + alpha**12*(t** &
      15/287332531200._qp + t**13/3150576000._qp) + alpha**11*(t**13/1260230400 &
      + t**11/111196800) + alpha**10*(t**13/1260230400 + t**11/18532800 &
      ) + alpha**9*(t**11/7413120 + t**9/926640) + alpha**8*(t**11/ &
      7413120 + t**9/154440) + alpha**7*(t**9/61776 + t**7/11880) + &
      alpha**6*(t**9/61776 + t**7/1980) + alpha**5*(t**7/792 + t**5/270 &
      ) + alpha**4*(t**7/792 + t**5/45) + alpha**3*(t**5/18 + t**3/15) &
      + alpha**2*(t**5/18 + 2*t**3/5) + alpha*t**3 + t**3

end function

real(qp) function S4(alpha, t)
implicit none
real(qp), intent(in) :: alpha
real(qp), intent(in) :: t

S4 = alpha**13*(t**16/670442572800._qp + t**14/27935107200._qp) + alpha**12*(t** &
      16/670442572800._qp + t**14/6207801600._qp + t**12/1470268800) + alpha** &
      11*(t**14/2660486400._qp + t**12/147026880) + alpha**10*(t**14/ &
      2660486400._qp + t**12/32672640 + t**10/10810800) + alpha**9*(t**12/ &
      14002560 + t**10/1081080) + alpha**8*(t**12/14002560 + t**10/ &
      240240 + t**8/120120) + alpha**7*(t**10/102960 + t**8/12012) + &
      alpha**6*(t**10/102960 + 3*t**8/8008 + t**6/2310) + alpha**5*(t** &
      8/1144 + t**6/231) + alpha**4*(t**8/1144 + 3*t**6/154 + t**4/105 &
      ) + alpha**3*(t**6/22 + 2*t**4/21) + alpha**2*(t**6/22 + 3*t**4/7 &
      ) + alpha*t**4 + t**4

end function

end module
