program ewald_convN

! Runs a benchmark with regards to the number of atoms N.

use types, only: dp
use ewald_sums, only: ewald_box, ewald_fft1, fred2fcart, min_distance
use md, only: positions_fcc
use utils, only: assert, init_random
implicit none

integer :: natom, ntypat

real(dp), parameter :: E_reference = -1.9853035145151909_dp
real(dp) :: stress(6)
real(dp), allocatable :: xred(:, :), forces(:, :), q(:)
real(dp) :: L, E_ewald, E_ewald_fft
real(dp) :: t1, t2, t3, t4
integer :: i, u

open(newunit=u, file="ewald_convN.txt", status="replace")

do i = 1, 6
    natom = 4*i**3
    ntypat = 1
    allocate(xred(3, natom), forces(3, natom), q(natom))
    call positions_fcc(xred, 1._dp)
    q = 1
    L = 2 / sqrt(3._dp)
    L = L * i
    print *, "half of xred min distance:", min_distance(xred*L, L)/2

    call cpu_time(t1)
    call ewald_box(L, xred*L, q, E_ewald, forces, stress)
    call cpu_time(t2)
    E_ewald = E_ewald / (natom/ntypat)

    call cpu_time(t3)
    call ewald_fft1(L, xred*L, q, 16*i, 0.40824_dp, E_ewald_fft)
    call cpu_time(t4)
    E_ewald_fft = E_ewald_fft / (natom/ntypat)
    print *, i, natom, t2-t1, t4-t3
    print *, E_ewald, E_ewald_fft
    print *, abs(E_ewald - E_reference), abs(E_ewald_fft-E_reference)
    write(u, *) i, natom, t2-t1, t4-t3, E_ewald, E_ewald_fft

    deallocate(xred, forces, q)
end do

close(u)

end program
