program test_gaussians_molecules

! [1] Szabo & Ostlund. Modern Quantum Chemistry (1996)

use types, only: dp
use drivers, only: rhf_gauss
use constants, only: ang2bohr
use mbpt, only: transform_int2, mbpt2, mbpt3
use utils, only: assert
use gf, only: find_poles, find_pole_diag
implicit none
real(dp), allocatable :: C(:, :), lam(:), int2(:), moint2(:)
real(dp) :: E0, E2, E3, Etot, Enuc, Ekin
! The difference between PyQuante C implementation and libint2 in the total
! energy for the molecules tested in this file should be below this eps:
real(dp), parameter :: eps = 1e-7_dp
integer :: Nscf, Nelec
real(dp) :: alpha, tolE, tolP, Egreen

tolE = 1e-10_dp
tolP = 1e-8_dp
alpha = 1._dp
Nscf = 40
! H2
call rhf_gauss([1, 1], reshape( &
    [0._dp, 0._dp, 0._dp, &
     0._dp, 0._dp, 1.4_dp], &
    [3, 2]), Nscf, alpha, tolE, tolP, C, lam, int2, E0, Enuc, Ekin, Etot)
allocate(moint2(size(int2)))
moint2 = transform_int2(int2, C)
print *, "Calculating MPBT 2"
Nelec = 2
E2 = mbpt2(moint2, lam, Nelec/2)
print *, "Calculating MPBT 3"
E3 = mbpt3(moint2, lam, Nelec/2)
call print_energies()
Egreen = find_pole_diag(1, moint2, lam, Nelec/2, 200, 1e-10_dp)
deallocate(moint2)
! PyQuante  -1.131 284
! [1]       -1.131     (Table 3.11, p. 192)
! GAMESS-UK -1.131 284 3467
call assert(abs(Etot - (-1.13128434_dp)) < eps)
! [1] Table 3.14, page 194:
call assert(abs(lam(1) - (-0.595_dp)) < 6e-4_dp)
! E2:
! [1]       -0.026 3         (Table 6.3, page 370)
! GAMESS-UK -0.026 341 79
call assert(abs(E2 - (-0.0263_dp)) < 6e-5_dp)
call assert(abs(E2+E3 - (-0.0319_dp)) < 6e-5_dp)
! E3:
! GAMESS-UK -0.00551598
call assert(abs(E3 - (-0.00551598_dp)) < 1e-8_dp)
! GF:
! [1]        0.598        (Table 7.2, page 406)
call assert(abs(Egreen - (-0.59832340_dp)) < 1e-8_dp)


! N2
call rhf_gauss([7, 7], reshape( &
    [0._dp, 0._dp, 0._dp, &
     0._dp, 0._dp, 2.074_dp], &
    [3, 2]), Nscf, alpha, tolE, tolP, C, lam, int2, E0, Enuc, Ekin, Etot)
call print_energies_hf()
! MPQC      -108.942 687
! [1]       -108.942         (Table 3.12, p. 192)
call assert(abs(Etot - (-108.94268654_dp)) < 1e-8_dp)
! [1] Table 3.16, page 196:
call assert(abs(lam(5) - (-0.630_dp)) < 6e-4_dp)
call assert(abs(lam(6) - (-0.612_dp)) < 6e-4_dp)
call assert(abs(lam(7) - (-0.612_dp)) < 6e-4_dp)

! CO
call rhf_gauss([8, 6], reshape( &
    [0._dp, 0._dp, 0._dp, &
     0._dp, 0._dp, 2.132_dp], &
    [3, 2]), Nscf, 0.8_dp, tolE, tolP, C, lam, int2, E0, Enuc, Ekin, Etot)
call print_energies_hf()
! [1]       -112.737         (Table 3.12, p. 192)
call assert(abs(Etot - (-112.73732121_dp)) < 1e-8_dp)
! [1] Table 3.15, page 195:
call assert(abs(lam(6) - (-0.633_dp)) < 6e-4_dp)
call assert(abs(lam(7) - (-0.548_dp)) < 6e-4_dp)

! CH4  Geometry 1
call rhf_gauss([6, 1, 1, 1, 1], reshape( &
    [0.000000_dp,     0.000000_dp,     0.000000_dp, &
    2.050000_dp,     0.000000_dp,     0.000000_dp, &
   -0.683292_dp,     1.932773_dp,     0.000000_dp, &
   -0.683292_dp,    -0.966387_dp,    -1.673831_dp, &
   -0.683292_dp,    -0.966387_dp,     1.673831_dp], &
    [3, 5]), Nscf, alpha, tolE, tolP, C, lam, int2, E0, Enuc, Ekin, Etot)
allocate(moint2(size(int2)))
moint2 = transform_int2(int2, C)
print *, "Calculating MPBT 2"
Nelec = 10
E2 = mbpt2(moint2, lam, Nelec/2)
print *, "Calculating MPBT 3"
E3 = mbpt3(moint2, lam, Nelec/2)
call print_energies()
Egreen = find_pole_diag(5, moint2, lam, Nelec/2, 200, 1e-10_dp)
deallocate(moint2)
! MPQC      -40.201 700
! [1]       -40.202         (Table 3.13, p. 192)
! GAMESS-UK -40.201 700 347 4
! PyQuante vs libint: 3.4e-8
call assert(abs(Etot - (-40.20170036_dp)) < 1e-8_dp)
! [1] Table 3.17, page 198 says -0.543, probably a typo
! GAMESS-UK: -0.54451009
call assert(abs(lam(5) - (-0.545_dp     )) < 6e-4_dp)
call assert(abs(lam(5) - (-0.54451008_dp)) < 1e-8_dp)
! E2:
! MPQC       -0.162 925  ! Wrong
! GAMESS-UK: -0.168 155 09
call assert(abs(E2 - (-0.16815509_dp)) < 1e-7_dp)
! GAMESS-UK: -0.018 192 43
call assert(abs(E3 - (-0.01819243_dp)) < 1e-7_dp)
! GF:
! [1]        0.510        (Table 7.4, page 407)
call assert(abs(Egreen - (-0.51384473_dp)) < 1e-8_dp)

! CH4  Geometry 2
call rhf_gauss([6, 1, 1, 1, 1], reshape( &
    [ 0.00000000_dp,     0.00000000_dp,     0.00000000_dp, &
    0.62558332_dp,    -0.62558332_dp,     0.62558332_dp, &
    -0.62558332_dp,     0.62558332_dp,     0.62558332_dp, &
    0.62558332_dp,     0.62558332_dp,    -0.62558332_dp, &
    -0.62558332_dp,    -0.62558332_dp,    -0.62558332_dp], &
    [3, 5]) * ang2bohr, Nscf, alpha, tolE, tolP, C, lam, int2, E0, Enuc, Ekin, Etot)
allocate(moint2(size(int2)))
moint2 = transform_int2(int2, C)
print *, "Calculating MPBT 2"
Nelec = 10
E2 = mbpt2(moint2, lam, Nelec/2)
call print_energies(skip3=.true.)
deallocate(moint2)
! GAMESS-UK: -40.201 704 77
call assert(abs(Etot - (-40.20170480_dp)) < 1e-8_dp)
! GAMESS-UK:  -0.168 150 19
call assert(abs(E2 - (-0.16814994_dp)) < 1e-6_dp)

! NH3
call rhf_gauss([7, 1, 1, 1], reshape( &
    [0.000000_dp,     0.000000_dp,     0.000000_dp, &
    1.913000_dp,     0.000000_dp,     0.000000_dp, &
   -0.548761_dp,     1.832602_dp,     0.000000_dp, &
   -0.548761_dp,    -0.916301_dp,    -1.587080_dp], &
    [3, 4]), Nscf, alpha, tolE, tolP, C, lam, int2, E0, Enuc, Ekin, Etot)
allocate(moint2(size(int2)))
moint2 = transform_int2(int2, C)
Nelec = 10
call print_energies()
Egreen = find_pole_diag(5, moint2, lam, Nelec/2, 200, 1e-10_dp)
deallocate(moint2)
! [1]       -56.195         (Table 3.13, p. 192)
! PyQuante vs libint: 5.0e-8
call assert(abs(Etot - (-56.19457246_dp)) < 1e-8_dp)
! [1] Table 3.17, page 198 says -0.421, but it is probably a typo.
call assert(abs(lam(5) - (-0.41491598_dp)) < 1e-8_dp)
! GF:
! [1]        0.353        (Table 7.4, page 407)
call assert(abs(Egreen - (-0.35357570_dp)) < 1e-8_dp)

! H2O
call rhf_gauss([8, 1, 1], reshape( &
    [0.000000_dp,     0.000000_dp,     0.000000_dp, &
    1.809000_dp,     0.000000_dp,     0.000000_dp, &
   -0.453549_dp,     1.751221_dp,     0.000000_dp], &
    [3, 3]), Nscf, alpha, tolE, tolP, C, lam, int2, E0, Enuc, Ekin, Etot)
allocate(moint2(size(int2)))
moint2 = transform_int2(int2, C)
print *, "Calculating MPBT 2"
Nelec = 10
E2 = mbpt2(moint2, lam, Nelec/2)
print *, "Calculating MPBT 3"
E3 = mbpt3(moint2, lam, Nelec/2)
call print_energies()
Egreen = find_pole_diag(5, moint2, lam, Nelec/2, 200, 1e-10_dp)
deallocate(moint2)
! MPQC:     -76.023 159
! [1]       -76.023         (Table 3.13, p. 192)
! GAMESS-UK -76.023 158 7398
! PyQuante vs libint: 8.0e-8
call assert(abs(Etot - (-76.02315869_dp)) < 1e-8_dp)
! GAMESS-UK: -0.49714239
! [1] Table 3.17, page 198:
call assert(abs(lam(5) - (-0.497_dp)) < 6e-4_dp)
! MPQC:      -0.196 587     ! Wrong
! GAMESS-UK: -0.199 259 95
call assert(abs(E2 - (-0.19925995_dp)) < 1e-6_dp)
! GAMESS-UK: -0.006 118 95
call assert(abs(E3 - (-0.00611895_dp)) < 1e-6_dp)
! GF:
! [1]        0.395        (Table 7.4, page 407), this is only one iteration
call assert(abs(Egreen - (-0.40499448_dp)) < 1e-8_dp)

! FH
call rhf_gauss([9, 1], reshape( &
    [0.000000_dp,     0.000000_dp,     0.000000_dp, &
    1.733_dp,     0.000000_dp,     0.000000_dp], &
    [3, 2]), Nscf, alpha, tolE, tolP, C, lam, int2, E0, Enuc, Ekin, Etot)
allocate(moint2(size(int2)))
moint2 = transform_int2(int2, C)
Nelec = 10
call print_energies()
Egreen = find_pole_diag(5, moint2, lam, Nelec/2, 200, 1e-10_dp)
deallocate(moint2)
! [1]       -100.011         (Table 3.13, p. 192)
! PyQuante vs libint: 9.0e-8
call assert(abs(Etot - (-100.01134814_dp)) < 1e-8_dp)
! [1] Table 3.17, page 198:
call assert(abs(lam(4) - (-0.627_dp)) < 6e-4_dp)
call assert(abs(lam(5) - (-0.627_dp)) < 6e-4_dp)
! GF:
! [1]        0.509        (Table 7.4, page 407), this is only one iteration
call assert(abs(Egreen - (-0.51894290_dp)) < 2e-8_dp)
! Note: Egreen with PyQuante integrals gives: -0.51894291
! with libint integrals, it gives           : -0.51894289

contains

subroutine print_energies_hf()
print *, "HF results:"
print "(a,es10.2)", "        EKIN+EHF (a.u.):", Ekin + Etot
print "(a,f18.8)", "  KINETIC ENERGY (a.u.):", Ekin
print "(a,f18.8)", "HF ATOMIC ENERGY (a.u.):", Etot
end subroutine

subroutine print_energies(skip3)
logical, optional, intent(in) :: skip3
call print_energies_hf()
print *
print *, "MBPT results:"
print "(' E0+E1 (HF)    = ',f15.8)", Etot
print "(' E2    (MBPT2) = ',f15.8)", E2
if (.not. present(skip3)) print "(' E3    (MBPT3) = ',f15.8)", E3
print "(' E0+E1+E2      = ',f15.8)", Etot + E2
if (.not. present(skip3)) print "(' E0+E1+E2+E3   = ',f15.8)", Etot + E2 + E3
end subroutine

end program
