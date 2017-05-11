module efermi
use types, only: dp
use optimize, only: bisect
implicit none
private
public fermi_dirac_smearing

contains

pure subroutine eval_fermi(eps, mu, sigma, occ)
real(dp), intent(in) :: eps(:), mu, sigma
real(dp), intent(out) :: occ(:)
occ = 1/(exp((eps-mu)/sigma) + 1)
end subroutine

subroutine fermi_dirac_smearing(eigs, sigma, Z, E_fermi, occ)
real(dp), intent(in) :: eigs(:), sigma, Z
real(dp), intent(out) :: E_fermi, occ(:)
integer :: nband
!E_fermi = 0.1708_dp
nband = size(eigs)
E_fermi = bisect(f, eigs(1), eigs(nband), 1e-12_dp)
call eval_fermi(eigs, E_fermi, sigma, occ)
occ = 2*occ ! double occupancy

contains

    real(dp) function f(x)
    real(dp), intent(in) :: x
    real(dp) :: occ(nband)
    call eval_fermi(eigs, x, sigma, occ)
    occ = 2*occ ! double occupancy
    f = Z - sum(occ)
    end function

end subroutine

end module

