module drivers
use types, only: dp
use scf, only: doscf, get_nuclear_energy
use basis, only: gaussints
use constants, only: Ha2eV
implicit none
private
public rhf_gauss

contains

subroutine rhf_gauss(atno, xyz, Nscf, alpha, tolE, tolP, C, lam, int2, &
        E0, Enuc, Etot)
integer, intent(in) :: atno(:) ! Atomic numbers (Z)
real(dp), intent(in) :: xyz(:, :) ! Atomic coordinates (3, :)
integer, intent(in) :: Nscf
real(dp), intent(in) :: alpha
real(dp), intent(in) :: tolE, tolP
real(dp), allocatable, intent(out) :: C(:, :), lam(:), int2(:)
real(dp), intent(out) :: E0, Enuc, Etot
real(dp), allocatable :: H(:, :), S(:, :), P(:, :), T(:, :), &
    V(:, :)
integer :: i

integer :: n

integer :: Nelec

call gaussints(atno, xyz, S, T, V, int2)
n = size(S, 1)
allocate(P(n, n), C(n, n), lam(n))

Nelec = sum(atno)

H = T + V

print *, "SCF cycle:"
call doscf(H, int2, S, Nelec/2, Nscf, tolE, tolP, alpha, C, P, lam, E0)
Enuc = get_nuclear_energy(atno, xyz)
Etot = E0 + Enuc
print *, "Orbital energies:"
print *, " i      E [Ha]      E [eV]"
do i = 1, size(lam)
    print "(i4,f12.6,f12.6)", i, lam(i), lam(i) * Ha2eV
end do
print *, "HF ENERGIES (a.u.):"
print *, "electronic =", E0
print *, "nuclear    =", Enuc
print *, "total      =", Etot
end subroutine

end module
