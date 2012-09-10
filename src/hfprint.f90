module hfprint
use types, only: dp
use constants, only: Ha2eV
implicit none
private
public printall, printlam, lmap

character :: lmap(0:5) = ["s", "p", "d", "f", "g", "h"]

contains

subroutine printall(nbfl, nl, zl, lam, C, Ekin, Etot)
integer, intent(in) :: nbfl(0:), nl(:, 0:)
real(dp), intent(in) :: zl(:, 0:), lam(:, 0:), C(:, :, 0:), Ekin, Etot
integer :: l, i
do l = 0, ubound(nbfl, 1)
    print *, lmap(l) // "-orbitals (in a.u.)"
    print *
    print "(a1, a14, 100(e14.6))", "N", "ZETA    ", lam(:nbfl(l), l)
    do i = 1, nbfl(l)
        print "(i1, 100(e14.6))", nl(i, l), zl(i, l), C(i, :nbfl(l), l)
    end do
    print *
end do
print "(a,e16.8)", "        EKIN+EHF (a.u.):", Ekin + Etot
print "(a,e16.8)", "  KINETIC ENERGY (a.u.):", Ekin
print "(a,e16.8)", "HF ATOMIC ENERGY (a.u.):", Etot
end subroutine

subroutine printlam(nbfl, lam, Ekin, Etot)
integer, intent(in) :: nbfl(0:)
real(dp), intent(in) :: lam(:, 0:), Ekin, Etot
integer :: l, i
print *, "Orbital energies:"
print *, "  n   l          E [a.u.]            E [eV]"
do l = 0, ubound(nbfl, 1)
    do i = 1, nbfl(l)
        print "(i4,i4,f20.8,f20.8)", i, l, lam(i, l), lam(i, l) * Ha2eV
    end do
end do
print "(a,es10.2)", "        EKIN+EHF (a.u.):", Ekin + Etot
print "(a,f18.8)", "  KINETIC ENERGY (a.u.):", Ekin
print "(a,f18.8)", "HF ATOMIC ENERGY (a.u.):", Etot
end subroutine

end module
