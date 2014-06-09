program test_ewald3

! This test compares lattice energy calculated from the Ewald summation
! against the energy calculated using the Madelung constant for various lattice
! constants L. The agreement is essentially to machine precision.

! This test tests ewald_box().


use types, only: dp
use constants, only: ang2bohr, kJmol2Ha
use ewald_sums, only: ewald_box, fred2fcart
use utils, only: assert, init_random
implicit none

integer :: natom, ntypat
! Various NaCl lattice constants in A
real(dp), parameter :: Llist(*) = [5.6402_dp, 5.5_dp, 4.5_dp, 6.5_dp, 10._dp]
! The component of the correct force for the given L
real(dp), parameter :: fcorrect(*) = [8.9613425608948444e-2_dp, &
    9.4240310569359220e-2_dp, 0.14077873554188225_dp, &
    6.7473831827766034e-2_dp, 2.8507693947231127e-2_dp]
real(dp), parameter :: stress1_correct(*) = [-0.62770548707869611_dp, &
    -0.64370627058568419_dp, -0.78675210849361421_dp, &
    -0.54467453664942533_dp, -0.35403844882212621_dp]
real(dp), parameter :: stress2_correct(*) = [-0.54226581445122812_dp, &
    -0.55608866303051230_dp, -0.67966392148173749_dp, &
    -0.47053656102581831_dp, -0.30584876466678168_dp]
real(dp), parameter :: stress3_correct(*) = [-0.43722866784288522_dp, &
    -0.44837402406680815_dp, -0.54801269608165426_dp, &
    -0.37939340497960677_dp, -0.24660571323674413_dp]

real(dp) :: stress(6)
real(dp), allocatable :: xred(:, :), forces(:, :), q(:)
real(dp) :: L, alpha, E_ewald, E_madelung
integer :: i, j

alpha = 1.74756459463318219064_dp ! Madelung constant for NaCl

! Conventional cell:
natom = 8
ntypat = 2
allocate(xred(3, natom), forces(3, natom), q(natom))
! Cl^-
xred(:, 1) = [0, 0, 0]
xred(:, 2) = [1, 1, 0] / 2._dp
xred(:, 3) = [1, 0, 1] / 2._dp
xred(:, 4) = [0, 1, 1] / 2._dp
! Na^+
xred(:, 5) = [1, 1, 1] / 2._dp
xred(:, 6) = [1, 0, 0] / 2._dp
xred(:, 7) = [0, 1, 0] / 2._dp
xred(:, 8) = [0, 0, 1] / 2._dp
q = [-1, -1, -1, -1, 1, 1, 1, 1]

do i = 1, size(Llist)
    L = Llist(i) * ang2bohr
    call ewald_box(L, xred*L, q, E_ewald, forces, stress)

    E_ewald = E_ewald / (natom/ntypat)
    E_madelung = -2*alpha/L

    print *, "a =", L/ang2bohr*100, "pm"
    print *, "Madelung:", E_madelung / kJmol2Ha, "kJ/mol"
    print *, "Ewald:   ", E_ewald / kJmol2Ha, "kJ/mol"
    print *, "error:   ", abs(E_ewald - E_madelung), "a.u."
    call assert(abs(E_ewald - E_madelung) < 1e-14_dp)
    call assert(all(abs(forces) < 1e-17_dp))
    call assert(all(abs(stress - [stress3_correct(i), stress3_correct(i), &
        stress3_correct(i), 0._dp, 0._dp, 0._dp]) < 1e-10_dp))
end do

print *, "--------"

! Madelung constant for NaCl, where the diagonal Na atom is moved by 3/8
! towards the Cl atom
alpha = 2.5088837163575413_dp

! Cl^-
xred(:, 1) = [0, 0, 0]
xred(:, 2) = [1, 1, 0] / 2._dp
xred(:, 3) = [1, 0, 1] / 2._dp
xred(:, 4) = [0, 1, 1] / 2._dp
! Na^+
xred(:, 5) = [1, 1, 1] / 2._dp
xred(:, 6) = [1, 0, 0] / 2._dp
xred(:, 7) = [0, 1, 0] / 2._dp
xred(:, 8) = [0, 0, 1] / 2._dp
xred(:, 5:8) = xred(:, 5:8) - spread([3, 3, 3] / 8._dp, 2, 4)
xred(:, 6:8) = xred(:, 6:8) - floor(xred(:, 6:8))

do i = 1, size(Llist)
    L = Llist(i) * ang2bohr
    call ewald_box(L, xred*L, q, E_ewald, forces, stress)

    E_ewald = E_ewald / (natom/ntypat)
    E_madelung = -2*alpha/L

    print *, "a =", L/ang2bohr*100, "pm"
    print *, "Madelung:", E_madelung / kJmol2Ha, "kJ/mol"
    print *, "Ewald:   ", E_ewald / kJmol2Ha, "kJ/mol"
    print *, "error:   ", abs(E_ewald - E_madelung), "a.u."
    call assert(abs(E_ewald - E_madelung) < 1e-14_dp)
    call assert(all(abs(stress - [stress1_correct(i), stress1_correct(i), &
        stress1_correct(i), stress2_correct(i), stress2_correct(i), &
        stress2_correct(i)]) < 1e-10_dp))
    do j = 1, 4
        call assert(all(abs(forces(:, j) - &
            [-fcorrect(i), -fcorrect(i), -fcorrect(i)]) < 1e-10_dp))
    end do
    do j = 5, 8
        call assert(all(abs(forces(:, j) - &
            [fcorrect(i), fcorrect(i), fcorrect(i)]) < 1e-10_dp))
    end do
end do
end program
