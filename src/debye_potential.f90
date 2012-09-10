module debye_potential
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

real(qp) function V0(D, r1, r2)
implicit none
real(qp), intent(in) :: D
real(qp), intent(in) :: r1
real(qp), intent(in) :: r2

V0 = (-D*exp(-(r1 + r2)/D) + D*exp(-(r1 - r2)/D))/(2*r1*r2)

end function

real(qp) function V1(D, r1, r2)
implicit none
real(qp), intent(in) :: D
real(qp), intent(in) :: r1
real(qp), intent(in) :: r2

V1 = 3*D*(-D**2*exp(2*r2/D) + D**2 - D*r1*exp(2*r2/D) + D*r1 + D*r2*exp( &
      2*r2/D) + D*r2 + r1*r2*exp(2*r2/D) + r1*r2)*exp((-r1 - r2)/D)/(2* &
      r1**2*r2**2)

end function

real(qp) function V2(D, r1, r2)
implicit none
real(qp), intent(in) :: D
real(qp), intent(in) :: r1
real(qp), intent(in) :: r2

V2 = 5*D*(9*D**4*exp(2*r2/D) - 9*D**4 + 9*D**3*r1*exp(2*r2/D) - 9*D**3* &
      r1 - 9*D**3*r2*exp(2*r2/D) - 9*D**3*r2 + 3*D**2*r1**2*exp(2*r2/D &
      ) - 3*D**2*r1**2 - 9*D**2*r1*r2*exp(2*r2/D) - 9*D**2*r1*r2 + 3*D &
      **2*r2**2*exp(2*r2/D) - 3*D**2*r2**2 - 3*D*r1**2*r2*exp(2*r2/D) - &
      3*D*r1**2*r2 + 3*D*r1*r2**2*exp(2*r2/D) - 3*D*r1*r2**2 + r1**2*r2 &
      **2*exp(2*r2/D) - r1**2*r2**2)*exp((-r1 - r2)/D)/(2*r1**3*r2**3)

end function

real(qp) function V3(D, r1, r2)
implicit none
real(qp), intent(in) :: D
real(qp), intent(in) :: r1
real(qp), intent(in) :: r2

V3 = 7*D*(-225*D**6*exp(2*r2/D) + 225*D**6 - 225*D**5*r1*exp(2*r2/D) + &
      225*D**5*r1 + 225*D**5*r2*exp(2*r2/D) + 225*D**5*r2 - 90*D**4*r1 &
      **2*exp(2*r2/D) + 90*D**4*r1**2 + 225*D**4*r1*r2*exp(2*r2/D) + &
      225*D**4*r1*r2 - 90*D**4*r2**2*exp(2*r2/D) + 90*D**4*r2**2 - 15*D &
      **3*r1**3*exp(2*r2/D) + 15*D**3*r1**3 + 90*D**3*r1**2*r2*exp(2*r2 &
      /D) + 90*D**3*r1**2*r2 - 90*D**3*r1*r2**2*exp(2*r2/D) + 90*D**3* &
      r1*r2**2 + 15*D**3*r2**3*exp(2*r2/D) + 15*D**3*r2**3 + 15*D**2*r1 &
      **3*r2*exp(2*r2/D) + 15*D**2*r1**3*r2 - 36*D**2*r1**2*r2**2*exp(2 &
      *r2/D) + 36*D**2*r1**2*r2**2 + 15*D**2*r1*r2**3*exp(2*r2/D) + 15* &
      D**2*r1*r2**3 - 6*D*r1**3*r2**2*exp(2*r2/D) + 6*D*r1**3*r2**2 + 6 &
      *D*r1**2*r2**3*exp(2*r2/D) + 6*D*r1**2*r2**3 + r1**3*r2**3*exp(2* &
      r2/D) + r1**3*r2**3)*exp((-r1 - r2)/D)/(2*r1**4*r2**4)

end function

real(qp) function V4(D, r1, r2)
implicit none
real(qp), intent(in) :: D
real(qp), intent(in) :: r1
real(qp), intent(in) :: r2

V4 = 9*D*(11025*D**8*exp(2*r2/D) - 11025*D**8 + 11025*D**7*r1*exp(2*r2/D &
      ) - 11025*D**7*r1 - 11025*D**7*r2*exp(2*r2/D) - 11025*D**7*r2 + &
      4725*D**6*r1**2*exp(2*r2/D) - 4725*D**6*r1**2 - 11025*D**6*r1*r2* &
      exp(2*r2/D) - 11025*D**6*r1*r2 + 4725*D**6*r2**2*exp(2*r2/D) - &
      4725*D**6*r2**2 + 1050*D**5*r1**3*exp(2*r2/D) - 1050*D**5*r1**3 - &
      4725*D**5*r1**2*r2*exp(2*r2/D) - 4725*D**5*r1**2*r2 + 4725*D**5* &
      r1*r2**2*exp(2*r2/D) - 4725*D**5*r1*r2**2 - 1050*D**5*r2**3*exp(2 &
      *r2/D) - 1050*D**5*r2**3 + 105*D**4*r1**4*exp(2*r2/D) - 105*D**4* &
      r1**4 - 1050*D**4*r1**3*r2*exp(2*r2/D) - 1050*D**4*r1**3*r2 + &
      2025*D**4*r1**2*r2**2*exp(2*r2/D) - 2025*D**4*r1**2*r2**2 - 1050* &
      D**4*r1*r2**3*exp(2*r2/D) - 1050*D**4*r1*r2**3 + 105*D**4*r2**4* &
      exp(2*r2/D) - 105*D**4*r2**4 - 105*D**3*r1**4*r2*exp(2*r2/D) - &
      105*D**3*r1**4*r2 + 450*D**3*r1**3*r2**2*exp(2*r2/D) - 450*D**3* &
      r1**3*r2**2 - 450*D**3*r1**2*r2**3*exp(2*r2/D) - 450*D**3*r1**2* &
      r2**3 + 105*D**3*r1*r2**4*exp(2*r2/D) - 105*D**3*r1*r2**4 + 45*D &
      **2*r1**4*r2**2*exp(2*r2/D) - 45*D**2*r1**4*r2**2 - 100*D**2*r1** &
      3*r2**3*exp(2*r2/D) - 100*D**2*r1**3*r2**3 + 45*D**2*r1**2*r2**4* &
      exp(2*r2/D) - 45*D**2*r1**2*r2**4 - 10*D*r1**4*r2**3*exp(2*r2/D) &
      - 10*D*r1**4*r2**3 + 10*D*r1**3*r2**4*exp(2*r2/D) - 10*D*r1**3*r2 &
      **4 + r1**4*r2**4*exp(2*r2/D) - r1**4*r2**4)*exp((-r1 - r2)/D)/(2 &
      *r1**5*r2**5)

end function

end module
