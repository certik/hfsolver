module constants
use types, only: dp
implicit none
private
public pi, e_, i_, bohr2ang, ang2bohr, Ha2eV

! Constants contain more digits than double precision, so that
! they are rounded correctly. Single letter constants contain underscore so
! that they do not clash with user variables ("e" and "i" are frequently used as
! loop variables)
real(dp), parameter :: pi    = 3.1415926535897932384626433832795_dp
real(dp), parameter :: e_    = 2.7182818284590452353602874713527_dp
complex(dp), parameter :: i_ = (0, 1)

real(dp), parameter :: Na = 6.02214129e23_dp ! Avogadro constant
! Standard uncertainty      0.00000027 (Source: 2010 CODATA)
real(dp), parameter :: Ha2eV = 27.21138505_dp    ! 1 Ha = (1 * Ha2eV) eV
! Standard uncertainty is       0.00000060 eV (Source: 2010 CODATA)

real(dp), parameter :: J2Ha = 2.29371248e17_dp  ! 1 J = (1 * J2Ha) Ha
! Standard uncertainty is     0.00000010 eV (Source: 2010 CODATA)

! Covert Ha to cm^{-1} (energy equivalent)
! 1 Ha = (1 * Ha2invcm) cm^{-1}
real(dp), parameter :: Ha2invcm = 219474.6313708_dp
! Standard uncertainty is              0.0000011 cm^{-1} (Source: 2010 CODATA)

! Conversion between Bohr (Hartree atomic units) and Angstrom
real(dp), parameter :: bohr2ang = 0.529177249_dp
real(dp), parameter :: ang2bohr = 1 / bohr2ang


! Converts eV to kcal/mol, it is assumed, that the number in eV is given per
! molecule.  1 eV = (1 * eV2kcalmol) kcal/mol
real(dp), parameter :: kcalmol2kJmol = 4.184_dp
real(dp), parameter :: kJmol2Ha = 1000 * J2Ha / Na

end module
