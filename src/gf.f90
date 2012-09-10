module gf
use types, only: dp
use scf, only: ijkl2intindex
use linalg, only: eigh, inv
use utils, only: str, stop_error
use constants, only: i_, pi
use quadrature, only: gauss_pts, gauss_wts
implicit none
private
public find_pole_diag, find_poles, greens_function_cplx, total_energy, &
    plot_poles

contains

real(dp) function find_pole_diag(i, moint2, lam, Noccupied, Nscf, tolE) &
        result(E)
! Solves the equation (7.44)
integer, intent(in) :: i, Noccupied, Nscf
real(dp), intent(in) :: moint2(:), lam(:), tolE
real(dp) :: Eold, dE
integer :: it
E = lam(i)
print *, "Green's function SCF:"
print "('    ',i4,'   E = ',f15.8)", 0, E
do it = 1, Nscf
    Eold = E
    E = lam(i) + self_energy(i, i, E, moint2, lam, Noccupied)
    dE = abs(E - Eold)
    print "('    ', i4, '   E = ', f15.8, '   dE = ', es10.2)", it, E, dE
    if (dE < tolE) return
end do
call stop_error("Green's function SCF did not converge in " // str(Nscf) // &
    " iterations.")
end function

subroutine plot_poles(moint2, lam, Noccupied)
integer, intent(in) :: Noccupied
real(dp), intent(in) :: moint2(:), lam(:)
real(dp) :: Eall(size(lam)), A(size(lam), size(lam)), S(size(lam), size(lam)), &
    C(size(lam), size(lam))
real(dp) :: E
integer :: i, j, u
S = 0
do i = 1, size(lam)
    S(i, i) = 1
end do
Eall = lam
E = 0
open(newunit=u, file="log.txt", status="replace")
do while (E > -30)
    do i = 1, size(lam)
        do j = 1, size(lam)
            A(i, j) = self_energy(i, j, E, moint2, lam, Noccupied)
        end do
    end do
    do i = 1, size(lam)
        A(i, i) = A(i, i) + lam(i)
    end do
    call eigh(A, S, Eall, C)
    print *, E, Eall
    write(u, *) E, Eall
    E = E - 0.1_dp
end do
close(u)
end subroutine

function find_poles(moint2, lam, Noccupied, Nscf) result(Eall)
! Solves the equation (7.42)
integer, intent(in) :: Noccupied, Nscf
real(dp), intent(in) :: moint2(:), lam(:)
real(dp) :: Eall(size(lam)), A(size(lam), size(lam)), S(size(lam), size(lam)), &
    C(size(lam), size(lam))
real(dp) :: E, Eold
integer :: it, i, j, tmp(1)
S = 0
do i = 1, size(lam)
    S(i, i) = 1
end do
Eall = lam
print *, "Green's function SCF:"
print *, 0, Eall(:10)
Eold = Eall(1) ! only the lowest for now
do it = 1, Nscf
    tmp = minloc((Eall - Eold)**2)
    i = tmp(1)
    E = Eall(i)
    do i = 1, size(lam)
        do j = 1, size(lam)
            A(i, j) = self_energy(i, j, E, moint2, lam, Noccupied)
        end do
    end do
    do i = 1, size(lam)
        A(i, i) = A(i, i) + lam(i)
    end do
    call eigh(A, S, Eall, C)
    !print *, it, Eall(:10)
    print *, it, E
    Eold = E
end do
!call stop_error("Green's function SCF did not converge in " // str(Nscf) // &
!    " iterations.")
end function

real(dp) function self_energy(i, j, E, moint2, lam, Noccupied) result(Eij)
! Implements equation (7.39)
integer, intent(in) :: i, j, Noccupied
real(dp), intent(in) :: E, moint2(:), lam(:)
integer :: a, b, r, s
real(dp) :: iras, jras, jsar, iabr, jabr, jbar
Eij = 0
do a = 1, Noccupied
    do r = Noccupied+1, size(lam)
        do s = Noccupied+1, size(lam)
            iras = moint2(ijkl2intindex(i, r, a, s))
            jras = moint2(ijkl2intindex(j, r, a, s))
            jsar = moint2(ijkl2intindex(j, s, a, r))
            Eij = Eij + iras*(2*jras-jsar)/(E + lam(a) - lam(r) - lam(s))
        end do
    end do
end do
do a = 1, Noccupied
    do b = 1, Noccupied
        do r = Noccupied+1, size(lam)
            iabr = moint2(ijkl2intindex(i, a, b, r))
            jabr = moint2(ijkl2intindex(j, a, b, r))
            jbar = moint2(ijkl2intindex(j, b, a, r))
            Eij = Eij + iabr*(2*jabr-jbar)/(E + lam(r) - lam(a) - lam(b))
        end do
    end do
end do
end function

complex(dp) function self_energy_cplx(i, j, E, moint2, lam, Noccupied) &
        result(Eij)
! Implements equation (7.39)
integer, intent(in) :: i, j, Noccupied
complex(dp), intent(in) :: E
real(dp), intent(in) :: moint2(:), lam(:)
integer :: a, b, r, s
real(dp) :: iras, jras, jsar, iabr, jabr, jbar
Eij = 0
do a = 1, Noccupied
    do r = Noccupied+1, size(lam)
        do s = Noccupied+1, size(lam)
            iras = moint2(ijkl2intindex(i, r, a, s))
            jras = moint2(ijkl2intindex(j, r, a, s))
            jsar = moint2(ijkl2intindex(j, s, a, r))
            Eij = Eij + iras*(2*jras-jsar)/(E + lam(a) - lam(r) - lam(s))
        end do
    end do
end do
do a = 1, Noccupied
    do b = 1, Noccupied
        do r = Noccupied+1, size(lam)
            iabr = moint2(ijkl2intindex(i, a, b, r))
            jabr = moint2(ijkl2intindex(j, a, b, r))
            jbar = moint2(ijkl2intindex(j, b, a, r))
            Eij = Eij + iabr*(2*jabr-jbar)/(E + lam(r) - lam(a) - lam(b))
        end do
    end do
end do
end function

function greens_function_cplx(E, moint2, lam, Noccupied) result(G)
integer, intent(in) :: Noccupied
real(dp), intent(in) :: moint2(:), lam(:)
complex(dp), intent(in) :: E
complex(dp) :: G(size(lam), size(lam))
integer :: i, j
do i = 1, size(lam)
    do j = 1, size(lam)
        G(i, j) = -self_energy_cplx(i, j, E, moint2, lam, Noccupied)
    end do
end do
! Uncomment this line to only get HF green's function (G0):
!G = 0
do i = 1, size(lam)
    G(i, i) = G(i, i) + E - lam(i)
end do
G = inv(G)
end function

function trace(A) result(t)
complex(dp), intent(in) :: A(:, :)
complex(dp) :: t
integer :: i
t = 0
do i = 1, size(A, 1)
    t = t + A(i, i)
end do
end function

real(dp) function total_energy(moint2, lam, Noccupied, Hcore) result(Etot_real)
integer, intent(in) :: Noccupied
real(dp), intent(in) :: moint2(:), lam(:)
real(dp), intent(in) :: Hcore(:, :)
complex(dp) :: Etot, Ntot
real(dp), parameter :: eps = 1e-12_dp
integer :: Nq
real(dp) :: w, h, e0
Etot = 0
w = 2
h = 0.5_dp
e0 = -1.76_dp
do Nq = 40, 64
    Ntot = 0
    ! Change Ef to Nf to calculate Ntot:
    Ntot = Ntot + integ(Ef, e0 - h*i_, e0 + h*i_, Nq)
    Ntot = Ntot + integ(Ef, e0 + h*i_, -w + h*i_, Nq)
    Ntot = Ntot + integ(Ef, -w + h*i_, -w - h*i_, Nq)
    Ntot = Ntot + integ(Ef, -w - h*i_, e0 - h*i_, Nq)
    Ntot = Ntot / (2*pi*i_)
    print *, "Nq =", Nq, "Ntot =", Ntot
end do
print *, "Etot = ", Etot

if (abs(aimag(Ntot)) < eps) print *, "Imaginary part of Ntot < 1e-12 (OK)"
if (abs(aimag(Etot)) < eps) print *, "Imaginary part of Etot < 1e-12 (OK)"

Etot_real = real(Etot, dp)
print *, "Etot (normalized) = ", Etot_real / real(Ntot, dp)


contains

    complex(dp) function Nf(x)
    complex(dp), intent(in) :: x
    Nf = 2*trace(greens_function_cplx(x, moint2, lam, Noccupied))
    end function

    complex(dp) function Ef(x)
    ! Calculates Tr((h1+z)*G(z))
    ! h1 is the sum of the one-electron kinetic and potential energy
    ! (i.e. electron-nuclei interaction), that is, in a basis:
    !     h1_basis = H^core = T + V = F - G.
    ! and over orbitals we simply multiply by coefficients:
    !     h1_{ij} = sum_{mu nu} C_{mu i} C_{nu j} * h1_basis_{mu nu}
    complex(dp), intent(in) :: x
    complex(dp) :: G(size(lam), size(lam))
    real(dp) :: h1(size(lam), size(lam))
    G = greens_function_cplx(x, moint2, lam, Noccupied)
    h1 = Hcore
    Ef = trace(matmul(h1, G)) + x * trace(G)
    end function

end function

complex(dp) function integ(f, a, b, Nq) result(s)
complex(dp), intent(in) :: a, b
interface
    complex(dp) function func(x)
    use types, only: dp
    implicit none
    complex(dp), intent(in) :: x
    end function
end interface
procedure(func) :: f
integer, intent(in) :: Nq
real(dp) :: xiq(Nq), wtq(Nq)
complex(dp) :: fq(Nq), jac, x
integer :: i
jac = (b-a)/2
xiq = gauss_pts(Nq)
wtq = gauss_wts(Nq)
do i = 1, Nq
    x = (xiq(i)+1) * jac + a
    fq(i) = f(x)
end do
s = sum(wtq * fq * jac)
end function


end module
