module rdirac_assembly
use types, only: dp
use sorting, only: sort
use constants, only: pi
use feutils, only: get_parent_nodes, get_parent_quad_pts_wts, phih, dphih, &
    get_nodes, define_connect, c2fullc => c2fullc2, get_quad_pts, fe2quad
use utils, only: assert, stop_error
use splines, only: iixmin, spline3pars, poly3, dpoly3
use mesh, only: meshexp
use linalg, only: eigh, eigvals
implicit none
private
public sfem

contains

subroutine sfem(Ne, p, Nq, L, Nb, kappa, a_, c, Z, eigs)
integer, intent(in) :: Ne, p, Nq
real(dp), intent(in) :: L, kappa, a_, c, Z
integer, intent(out) :: Nb
real(dp), allocatable, intent(out) :: eigs(:)

integer :: Nn
! xe(i) is the 'x' coordinate of the i-th mesh node
real(dp), allocatable :: xe(:)
real(dp), allocatable :: xin(:), xiq(:), wtq(:), A(:, :), B(:, :), sol(:, :), &
    phihq(:, :), dphihq(:, :), Vq(:,:), xn(:), &
    fullc(:), enrq(:,:,:), denrq(:,:,:), phipuq(:,:), dphipuq(:,:), xinpu(:), &
    xq(:,:)
integer, allocatable :: ib(:, :), in(:, :), ibenr(:,:,:)
real(dp) :: rc, En
integer :: i, j, iqx, u, Nenr, emin, emax, l_, relat

Nn = Ne*p+1
allocate(xe(Ne+1))
xe = meshexp(0._dp, L, a_, Ne) ! uniform mesh on [0, L]

allocate(xin(p+1), Vq(Nq,Ne), xq(Nq,Ne))
allocate(xinpu(2)) ! linear functions for PU
call get_parent_nodes(2, p, xin)
call get_parent_nodes(2, size(xinpu)-1, xinpu)
allocate(xiq(Nq), wtq(Nq))
call get_parent_quad_pts_wts(1, Nq, xiq, wtq)
allocate(xn(Nn))
call get_nodes(xe, xin, xn)
call get_quad_pts(xe, xiq, xq)
allocate(phihq(size(xiq), size(xin)))
allocate(phipuq(size(xiq), size(xinpu)))
allocate(dphipuq(size(xiq), size(xinpu)))
allocate(dphihq(size(xiq), size(xin)))
! Tabulate parent basis at quadrature points
forall(i=1:size(xiq), j=1:size(xin))  phihq(i, j) =  phih(xin, j, xiq(i))
forall(i=1:size(xiq), j=1:size(xin)) dphihq(i, j) = dphih(xin, j, xiq(i))
forall(i=1:size(xiq), j=1:size(xinpu))  phipuq(i, j) =  phih(xinpu, j, xiq(i))
forall(i=1:size(xiq), j=1:size(xinpu)) dphipuq(i, j) = dphih(xinpu, j, xiq(i))

allocate(in(p+1,Ne),ib(p+1,Ne))
call define_connect(1,1,Ne,p,in,ib)

Nb = maxval(ib)

allocate(A(Nb, Nb), B(Nb, Nb), sol(Nb, Nb), eigs(Nb))
allocate(fullc(Nn))

Vq = -Z/xq

call assemble_1d(xin, xe, ib, xiq, wtq, phihq, dphihq, Vq, A, B)
call eigh(A, B, eigs, sol)

print *, "n, energy, exact energy, error"
do i = 1, Nb
    if (kappa > 0) then
        l_ = kappa
        relat = 3
    else
        l_ = -kappa-1
        relat = 2
    end if
    En = E_nl(c, l_+i, l_, real(Z, dp), relat)
    eigs(i) = sqrt(eigs(i)) - c**2
    print "(i4, f30.8, f18.8, es12.2)", i, eigs(i), En, abs(eigs(i)-En)
end do
end subroutine

real(dp) function E_nl(c, n, l, Z, relat)
! Calculates exact energy for the radial Schroedinger/Dirac equations
real(dp), intent(in) :: c, Z ! speed of light in atomic units
integer, intent(in) :: n, l, relat
! quantum numbers (n, l), atomic number (z)
! relat == 0 ... Schroedinger equation
! relat == 2 ... Dirac equation, spin up
! relat == 3 ... Dirac equation, spin down

integer :: skappa
real(dp) :: beta
if (.not. (l >= 0)) call stop_error("'l' must be positive or zero")
if (.not. (n > l)) call stop_error("'n' must be greater than 'l'")
if (l == 0 .and. relat == 3) call stop_error("Spin must be up for l==0.")
if (relat == 0) then
    E_nl = - Z**2 / (2.0_dp * n**2)
else
    if (relat == 2) then
        skappa = -l - 1
    else
        skappa = -l
    end if
    beta = sqrt(skappa**2 - Z**2 / c**2)
    E_nl = c**2/sqrt(1+Z**2/(n + skappa + beta)**2/c**2) - c**2
end if
end function

subroutine assemble_1d(xin, nodes, ib, xiq, wtq, phihq, dphihq, Vq, Am, Bm)
! Assemble on a 1D uniform mesh
real(dp), intent(in):: xin(:), nodes(:), xiq(:), wtq(:), &
    phihq(:, :), dphihq(:, :), Vq(:,:)
integer, intent(in):: ib(:, :)
real(dp), intent(out):: Am(:,:), Bm(:, :)
real(dp), dimension(size(xiq), &
    size(xin)) :: phi_v, phi_dx
real(dp), dimension(size(xiq)) :: x, xp
integer :: Ne, p, e, i, j, iqx
real(dp) :: lx
integer :: ax, bx
real(dp) :: jacx, jac_det

Ne = size(nodes)-1
p = size(xin) - 1
! 1D shape functions
do ax = 1, p+1
    do iqx = 1, size(xiq)
        phi_v (iqx, ax) =  phihq(iqx, ax)
        phi_dx(iqx, ax) = dphihq(iqx, ax)
    end do
end do
Am=0; Bm=0
! Precalculate as much as possible:
lx = nodes(2) - nodes(1) ! Element size
jacx = lx/2
jac_det = abs(jacx)
xp = (xiq + 1) * jacx
phi_dx = phi_dx / jacx
do e = 1, Ne
    x = xp + nodes(e)
    do bx = 1, p+1
        j = ib(bx, e)
        if (j==0) cycle
        do ax = 1, p+1
            i = ib(ax, e)
            if (i == 0) cycle
            if (j > i) cycle
            Am(i,j) = Am(i,j) + sum(phi_dx(:, ax)*phi_dx(:, bx) &
                * jac_det * wtq) / 2
            Am(i,j) = Am(i,j) + sum(Vq(:,e) * &
                phi_v(:, ax)*phi_v(:, bx) &
                * jac_det * wtq)
            Bm(i,j) = Bm(i,j) + sum(( &
                phi_v(:, ax) * phi_v(:, bx) &
                * jac_det * wtq))
        end do
    end do
end do
do j = 1, size(Am, 2)
    do i = 1, j-1
        Am(i, j) = Am(j, i)
        Bm(i, j) = Bm(j, i)
    end do
end do
end subroutine


end module


!------------------------------------------------------------------------------

program rdirac

use types, only: dp
use rdirac_assembly, only: sfem
implicit none

integer :: Ne, p, Nq, DOFs, i, u
real(dp), allocatable :: eigs(:)
real(dp) :: L, a, c, kappa, Z

Ne = 4
p = 10
Nq = 64
L = 10
a = 1e4
Z = 92
kappa = 2
c = 137.03599907_dp
call sfem(Ne, p, Nq, L, DOFs, kappa, a, c, Z, eigs)
print *, "Ne:", Ne
print *, "p:", p
print *, "Nq:", Nq
print *, "DOFs:", DOFs
do i = 1, 6
    print *, i, eigs(i)
end do


end program
