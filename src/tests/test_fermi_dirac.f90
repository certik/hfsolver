program test_fermi_dirac
use types, only: dp
use utils, only: assert, loadtxt
use efermi, only: fermi_dirac_smearing
implicit none
real(dp), parameter :: Z = 32 ! total charge
real(dp), parameter :: E_fermi_ref = 0.1708_dp ! eV
real(dp), parameter :: sigma = 2.69_dp ! eV
real(dp), allocatable :: D(:,:), eigs(:), occ_ref(:), occ(:)
real(dp) :: E_fermi
integer :: nband, i
call loadtxt("data/occ1.dat", D)
nband = size(D,1)
allocate(eigs(nband), occ_ref(nband), occ(nband))
eigs = D(:, 2) ! eV
occ_ref = D(:, 3)
call fermi_dirac_smearing(eigs, sigma, Z, E_fermi, occ)
do i = nband, 1, -1
    print "(i4,f10.4,f10.5,f10.5)", i, eigs(i), occ_ref(i), occ(i)
end do
print *, "Charges (occ and occ_ref):"
print *, sum(occ), sum(occ_ref)
print *, "Error:", abs(sum(occ)-Z)
call assert(abs(sum(occ)-Z) < 1e-12_dp)
print *
print *, "Fermi energy (E_fermi and E_fermi_ref):", E_fermi, E_fermi_ref
print *, "Error:", abs(E_fermi - E_fermi_ref)
call assert(abs(E_fermi - E_fermi_ref) < 1e-4_dp)
end program
