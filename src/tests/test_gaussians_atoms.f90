program test_gaussians_atoms

! The paper [2] uses 6-31G* basis (we use 6-31G**):

! [2] Johnson, B. G., Gill, P. M. W., & Pople, J. A. (1993). The
! performance of a family of density functional methods. Journal of Chemical
! Physics, 98(7), 5612-5626.

use types, only: dp
use drivers, only: rhf_gauss
use mbpt, only: transform_int2, mbpt2, mbpt3
use gf, only: plot_poles
use utils, only: assert
implicit none
real(dp), allocatable :: C(:, :), lam(:), int2(:), moint2(:)
real(dp) :: E0, E2, E3, Etot, Enuc, Ekin
integer :: Nscf, Nelec
real(dp) :: alpha, tolE, tolP

tolE = 1e-10_dp
tolP = 1e-8_dp
alpha = 1._dp
Nscf = 40
! He
call rhf_gauss([2], reshape( &
    [0._dp, 0._dp, 0._dp], &
    [3, 1]), Nscf, alpha, tolE, tolP, C, lam, int2, E0, Enuc, Ekin, Etot)
allocate(moint2(size(int2)))
moint2 = transform_int2(int2, C)
print *, "Calculating MPBT 2"
Nelec = 2
E2 = mbpt2(moint2, lam, Nelec/2)
print *, "Calculating MPBT 3"
E3 = mbpt3(moint2, lam, Nelec/2)
call print_energies()
deallocate(moint2)
! MPQC:      -2.855 160
! Jaguar:    -2.855 160
! [2]:       -2.855 2
! GAMESS-UK: -2.855 160 4262
call assert(abs(Etot - (-2.8551604262_dp)) < 1e-10_dp)
! MPQC:      -0.025 477
! Jaguar:    -0.025 477
! [2]:       -0.011 2
! GAMESS-UK: -0.025 477 03
call assert(abs(E2 - (-0.02547703_dp)) < 1e-6_dp)
! GAMESS-UK: -0.005 442 24
call assert(abs(E3 - (-0.00544224_dp)) < 1e-6_dp)

! Be
call rhf_gauss([4], reshape( &
    [0._dp, 0._dp, 0._dp], &
    [3, 1]), Nscf, alpha, tolE, tolP, C, lam, int2, E0, Enuc, Ekin, Etot)
allocate(moint2(size(int2)))
moint2 = transform_int2(int2, C)
print *, "Calculating MPBT 2"
Nelec = 4
E2 = mbpt2(moint2, lam, Nelec/2)
print *, "Calculating MPBT 3"
E3 = mbpt3(moint2, lam, Nelec/2)
call print_energies()
! This is how to plot poles:
!call plot_poles(moint2, lam, Nelec/2)
deallocate(moint2)
! MPQC:      -14.272 316 ! Wrong
! [2]:       -14.566 9
! GAMESS-UK: -14.569 853 664 5
call assert(abs(Etot - (-14.56694436_dp)) < 1e-8_dp)
! MPQC:      -0.042 703 ! Wrong
! [2]:       -0.029 5
! GAMESS-UK: -0.027 205 00
call assert(abs(E2 - (-0.02950068_dp)) < 1e-6_dp)
! GAMESS-UK: -0.009 907 62
call assert(abs(E3 - (-0.01056728_dp)) < 1e-6_dp)

! Ne
call rhf_gauss([10], reshape( &
    [0._dp, 0._dp, 0._dp], &
    [3, 1]), Nscf, alpha, tolE, tolP, C, lam, int2, E0, Enuc, Ekin, Etot)
allocate(moint2(size(int2)))
moint2 = transform_int2(int2, C)
print *, "Calculating MPBT 2"
Nelec = 10
E2 = mbpt2(moint2, lam, Nelec/2)
print *, "Calculating MPBT 3"
E3 = mbpt3(moint2, lam, Nelec/2)
call print_energies()
deallocate(moint2)
! PyQuante vs libint: 9.7e-8
! MPQC:      -128.474 407
! GAMESS-UK: -128.474 406 519 9
! [2]:       -128.474 4
call assert(abs(Etot - (-128.4744065196_dp)) < 1e-10_dp)
! MPQC:        -0.150 316  ! Wrong
! [2]:         -0.151 8
! GAMESS-UK    -0.15176955
call assert(abs(E2 - (-0.15176955_dp)) < 1e-6_dp)
! GAMESS-UK     0.00015112
call assert(abs(E3 - (0.00015112_dp)) < 1e-6_dp)

! Mg
call rhf_gauss([12], reshape( &
    [0._dp, 0._dp, 0._dp], &
    [3, 1]), Nscf, alpha, tolE, tolP, C, lam, int2, E0, Enuc, Ekin, Etot)
allocate(moint2(size(int2)))
moint2 = transform_int2(int2, C)
print *, "Calculating MPBT 2"
Nelec = 12
E2 = mbpt2(moint2, lam, Nelec/2)
print *, "Calculating MPBT 3"
E3 = mbpt3(moint2, lam, Nelec/2)
call print_energies()
deallocate(moint2)
! PyQuante vs libint: 1.5e-7
! MPQC       -199.595 610 920 833
! GAMESS-UK  -199.595 610 917 0
call assert(abs(Etot - (-199.59561092_dp)) < 1e-8_dp)
! MPQC       -0.021 986 440 863
! GAMESS-UK  -0.028 462 86
call assert(abs(E2 - (-0.02846287_dp)) < 1e-6_dp)
! GAMESS-UK  -0.007 145 35
call assert(abs(E3 - (-0.00714535_dp)) < 1e-6_dp)

! Ar
call rhf_gauss([18], reshape( &
    [0._dp, 0._dp, 0._dp], &
    [3, 1]), Nscf, alpha, tolE, tolP, C, lam, int2, E0, Enuc, Ekin, Etot)
allocate(moint2(size(int2)))
moint2 = transform_int2(int2, C)
print *, "Calculating MPBT 2"
Nelec = 12
E2 = mbpt2(moint2, lam, Nelec/2)
print *, "Calculating MPBT 3"
E3 = mbpt3(moint2, lam, Nelec/2)
call print_energies()
deallocate(moint2)
! PyQuante vs libint: 4.0e-7
! MPQC:      -526.773 744 920 946
! GAMESS-UK: -526.773 744 921 0
call assert(abs(Etot - (-526.77374492_dp)) < 1e-8_dp)
! MPQC:        -0.137308007079
! GAMESS-UK:   -0.14624537
call assert(abs(E2 - (-0.06205272_dp)) < 1e-6_dp)
! GAMESS-UK:   -0.01112785
call assert(abs(E3 - (-0.02128214_dp)) < 1e-6_dp)

contains

subroutine print_energies(skip3)
logical, optional, intent(in) :: skip3
print *, "HF results:"
print "(a,es10.2)", "        EKIN+EHF (a.u.):", Ekin + Etot
print "(a,f18.8)", "  KINETIC ENERGY (a.u.):", Ekin
print "(a,f18.8)", "HF ATOMIC ENERGY (a.u.):", Etot
print *
print *, "MBPT results:"
print "(' E0+E1 (HF)    = ',f15.8)", Etot
print "(' E2    (MBPT2) = ',f15.8)", E2
if (.not. present(skip3)) print "(' E3    (MBPT3) = ',f15.8)", E3
print "(' E0+E1+E2      = ',f15.8)", Etot + E2
if (.not. present(skip3)) print "(' E0+E1+E2+E3   = ',f15.8)", Etot + E2 + E3
end subroutine

end program
