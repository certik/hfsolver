module efermi
use types, only: dp
use optimize, only: bisect
use utils, only: stop_error
implicit none
private
public fermi_dirac_smearing

contains

pure subroutine eval_fermi(eps, mu, sigma, occ)
real(dp), intent(in) :: eps(:), mu, sigma
real(dp), intent(out) :: occ(:)
real(dp) :: arg(size(eps))
arg = (eps-mu)/sigma
where (arg < 50)
    occ = 1/(exp(arg) + 1)
else where
    occ = 0
end where
end subroutine

subroutine fermi_dirac_smearing(eigs, sigma, Z, E_fermi, occ)
real(dp), intent(in) :: eigs(:), sigma, Z
real(dp), intent(out) :: E_fermi, occ(:)
integer :: nband, nband_min
nband = size(eigs)
nband_min = ceiling(Z/2) + 1
if (nband < nband_min) then
    print *, "Increase nband at least to nband =", nband_min
    call stop_error("nband too small.")
end if
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

