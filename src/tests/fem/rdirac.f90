module rdirac_assembly
use types, only: dp
use sorting, only: sort
use constants, only: pi
use feutils, only: get_parent_nodes, get_parent_quad_pts_wts, phih, dphih, &
    get_nodes, define_connect_n, c2fullc => c2fullc2, get_quad_pts, fe2quad
use utils, only: assert, stop_error
use splines, only: iixmin, spline3pars, poly3, dpoly3
use mesh, only: meshexp
use linalg, only: eigh, eigvals
implicit none
private
public sfem

contains

subroutine sfem(Ne, p, Nq, L, Nb, kappa, a_, c, Z, eigs, squared)
integer, intent(in) :: Ne, p, Nq, kappa
real(dp), intent(in) :: L, a_, c, Z
integer, intent(out) :: Nb
real(dp), allocatable, intent(out) :: eigs(:)
logical, intent(in) :: squared

integer :: Nn
! xe(i) is the 'x' coordinate of the i-th mesh node
real(dp), allocatable :: xe(:)
real(dp), allocatable :: xin(:), xiq(:), wtq(:), A(:, :), B(:, :), sol(:, :), &
    phihq(:, :), dphihq(:, :), Vq(:,:), xn(:), &
    fullc(:), phipuq(:,:), dphipuq(:,:), xinpu(:), &
    xq(:,:)
integer, allocatable :: ib(:,:,:), in(:,:,:)
real(dp) :: En
integer :: i, j, l_, relat

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

allocate(in(p+1,Ne,2),ib(p+1,Ne,2))
call define_connect_n([1,2],[1,2],Ne,p,2,in,ib)

Nb = maxval(ib)

allocate(A(Nb, Nb), B(Nb, Nb), sol(Nb, Nb), eigs(Nb))
allocate(fullc(Nn))

Vq = -Z/xq

if (squared) then
    call assemble_rdirac_squared(Vq,c,kappa,xin,xe,ib,xiq,wtq,A,B)
else
    call assemble_rdirac(Vq,c,kappa,xin,xe,ib,xiq,wtq,A,B)
end if
call eigh(A, B, eigs, sol)
j = 1
do while (eigs(j) < 0)
  j = j + 1
end do

print *, "n, energy, exact energy, error"
do i = 1, min(20,Nb)
    if (kappa > 0) then
        l_ = kappa
        relat = 3
    else
        l_ = -kappa-1
        relat = 2
    end if
    En = E_nl(c, l_+i, l_, real(Z, dp), relat)
    if (squared) then
        eigs(i) = sqrt(eigs(i+j-1)) - c**2
    else
        eigs(i) = eigs(i+j-1) - c**2
    end if
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

subroutine assemble_rdirac(V,c,kappa,xin,xe,ib,xiq,wtq,Am,Bm)
! forms system equation matrices corresponding to the Dirac eigenproblem
! subject to boundary conditions consistent with basis specified by ib
! Both Dirac components are on the same mesh.
real(dp), intent(in) :: V(:,:)   ! V(x) at quadrature points:
   ! V(i,j) = value at ith point in jth element
real(dp), intent(in) :: c ! speed of light
integer, intent(in) :: kappa     ! kappa
real(dp), intent(in):: xin(:)       ! parent basis nodes
real(dp), intent(in):: xe(:)
integer, intent(in):: ib(:,:,:)     ! basis connectivity: ib(i,j,k) = index of
    ! basis function associated with local basis function i of element j of
    ! mesh k. 0 = no associated basis fn.
real(dp), intent(in):: xiq(:)       ! quadrature points
real(dp), intent(in):: wtq(:)       ! quadrature weights
real(dp), intent(out):: Am(:,:)     ! LHS matrix: Am c = lam Bm c
real(dp), intent(out):: Bm(:,:)     ! RHS matrix: Am c = lam Bm c
integer Ne, Nb              ! number of elements, nodes, basis functions
integer p                     ! order of FE/SE basis
integer e                     ! element index
integer i,j                   ! basis fn indices
integer :: m1, m2
integer al,be                 ! "alpha", "beta": local basis fn indices
integer iq                    ! quadrature point index
real(dp) xa,xb                ! element boundary node coordinates
real(dp) jac                  ! Jacobian of transformation from parent coords xi in [-1,1]
                              ! to coords x in [xa,xb]: x = (xb-xa)/2 xi + (xb+xa)/2
real(dp) phihq(size(xiq),size(xin)),dphihq(size(xiq),size(xin))   ! parent basis fns and derivs
   ! at quadrature points: phihq(i,j) = value of jth function at ith quadrature point
real(dp) Vq(size(xiq))
real(dp) Vqx(size(xiq))
real(dp) hq(size(xiq))        ! integrand at quadrature points
real(dp) x(size(xiq))                    ! point in domain
real(dp) :: n, k
real(dp), dimension(size(xiq)) :: Bi, Bj, Bip, Bjp

n = 0
k = 1000

! initializations
Ne= size(xe) - 1
Nb=maxval(ib)
p=size(xin)-1
if (size(Am,1)/=Nb .or. size(Bm,1)/=Nb) call stop_error("Error: size of Am and/or Bm inconsistent with Nb.")
! tabulate parent basis and derivatives at quadrature points
do al = 1, p+1
   do iq =1, size(xiq)
      phihq(iq, al) = phih(xin, al, xiq(iq))
      dphihq(iq, al) = dphih(xin, al, xiq(iq))
   end do
end do

! accumulate Am and Bm matrix elements
Am=0; Bm=0
do m1 = 1, 2
    do m2 = 1, 2
        do e = 1, Ne
            xa = xe(e)
            xb = xe(e+1)
            jac = (xb-xa)/2  ! affine mapping
            x = (xiq(:)+1) * jac + xa ! affine mapping
            Vq = V(:, e)
            Vqx = Vq*x ! Potential times r
            ! compute matrix elements (integrals transformed to [-1,1])
            do be = 1, p+1
                j = ib(be, e, m2)
                if (j == 0) cycle    ! omit boundary basis fns for Dirichlet BCs
                do al = 1, p+1
                    i = ib(al, e, m1)
                    if (i == 0) cycle! omit boundary basis fns for Dirichlet BCs
                    Bi = phihq(:,al)
                    Bj = phihq(:,be)
                    Bip = dphihq(:,al)/jac
                    Bjp = dphihq(:,be)/jac
                    if (m1 == 1 .and. m2 == 1) then
                        hq = Bi*Bj*(Vq+c**2)
                    else if (m1 == 1 .and. m2 == 2) then
                        hq = Bi*c*(-Bjp + kappa/x*Bj)
                    else if (m1 == 2 .and. m2 == 1) then
                        hq = Bi*c*(+Bjp + kappa/x*Bj)
                    else if (m1 == 2 .and. m2 == 2) then
                        hq = Bi*Bj*(Vq-c**2)
                    end if
                    Am(i,j) = Am(i,j) + sum(wtq*hq*jac)
                    if (m1 == m2) then
                        Bm(i,j) = Bm(i,j) + sum(wtq*Bi*Bj*jac)
                    end if
                end do
            end do
        end do
    end do
end do
! check symmetry
print *, "Checking symmetry"
do j=1,Nb
    do i=1,j-1
        if (abs(Am(i,j)-Am(j,i)) / (max(Am(i,j), Am(j,i))+tiny(1._dp)) &
                > 1e-8_dp .and. abs(Am(i,j)-Am(j,i)) > 1e-8_dp) then
            print *, i, j, Am(i,j)-Am(j,i), Am(i,j), Am(j,i)
            call stop_error("Am not symmetric")
        end if
        if (abs(Bm(i,j)-Bm(j,i)) > 1e-12_dp) call stop_error("Bm not symmetric")
   end do
end do
end subroutine


subroutine assemble_rdirac_squared(V,c,kappa,xin,xe,ib,xiq,wtq,Am,Bm)
! forms system equation matrices corresponding to the Dirac eigenproblem
! subject to boundary conditions consistent with basis specified by ib
! Both Dirac components are on the same mesh.
real(dp), intent(in) :: V(:,:)   ! V(x) at quadrature points:
   ! V(i,j) = value at ith point in jth element
real(dp), intent(in) :: c ! speed of light
integer, intent(in) :: kappa     ! kappa
real(dp), intent(in):: xin(:)       ! parent basis nodes
real(dp), intent(in):: xe(:)
integer, intent(in):: ib(:,:,:)     ! basis connectivity: ib(i,j,k) = index of
    ! basis function associated with local basis function i of element j of
    ! mesh k. 0 = no associated basis fn.
real(dp), intent(in):: xiq(:)       ! quadrature points
real(dp), intent(in):: wtq(:)       ! quadrature weights
real(dp), intent(out):: Am(:,:)     ! LHS matrix: Am c = lam Bm c
real(dp), intent(out):: Bm(:,:)     ! RHS matrix: Am c = lam Bm c
integer Ne, Nb              ! number of elements, nodes, basis functions
integer p                     ! order of FE/SE basis
integer e                     ! element index
integer i,j                   ! basis fn indices
integer :: m1, m2
integer al,be                 ! "alpha", "beta": local basis fn indices
integer iq                    ! quadrature point index
real(dp) xa,xb                ! element boundary node coordinates
real(dp) jac                  ! Jacobian of transformation from parent coords xi in [-1,1]
                              ! to coords x in [xa,xb]: x = (xb-xa)/2 xi + (xb+xa)/2
real(dp) phihq(size(xiq),size(xin)),dphihq(size(xiq),size(xin))   ! parent basis fns and derivs
   ! at quadrature points: phihq(i,j) = value of jth function at ith quadrature point
real(dp) Vq(size(xiq))
real(dp) Vqx(size(xiq))
real(dp) hq(size(xiq))        ! integrand at quadrature points
real(dp) x(size(xiq))                    ! point in domain
real(dp) :: n, k
real(dp), dimension(size(xiq)) :: Bi, Bj, Bip, Bjp

n = 0
k = 1000

! initializations
Ne= size(xe) - 1
Nb=maxval(ib)
p=size(xin)-1
if (size(Am,1)/=Nb .or. size(Bm,1)/=Nb) call stop_error("Error: size of Am and/or Bm inconsistent with Nb.")
! tabulate parent basis and derivatives at quadrature points
do al = 1, p+1
   do iq =1, size(xiq)
      phihq(iq, al) = phih(xin, al, xiq(iq))
      dphihq(iq, al) = dphih(xin, al, xiq(iq))
   end do
end do

! accumulate Am and Bm matrix elements
Am=0; Bm=0
do m1 = 1, 2
    do m2 = 1, 2
        do e = 1, Ne
            xa = xe(e)
            xb = xe(e+1)
            jac = (xb-xa)/2  ! affine mapping
            x = (xiq(:)+1) * jac + xa ! affine mapping
            Vq = V(:, e)
            Vqx = Vq*x ! Potential times r
            ! compute matrix elements (integrals transformed to [-1,1])
            do be = 1, p+1
                j = ib(be, e, m2)
                if (j == 0) cycle    ! omit boundary basis fns for Dirichlet BCs
                do al = 1, p+1
                    i = ib(al, e, m1)
                    if (i == 0) cycle! omit boundary basis fns for Dirichlet BCs
                    Bi = phihq(:,al)
                    Bj = phihq(:,be)
                    Bip = dphihq(:,al)/jac
                    Bjp = dphihq(:,be)/jac
                    if (m1 == 1 .and. m2 == 1) then
                        hq = c**2*Bip*Bjp &
                            + Bi*Bj*((Vq+c**2)**2+c**2*(kappa*(kappa+1)/x**2))
                    else if (m1 == 1 .and. m2 == 2) then
                        hq = c*Vq*(+Bip*Bj-Bi*Bjp + 2*kappa/x*Bi*Bj)
                    else if (m1 == 2 .and. m2 == 1) then
                        hq = c*Vq*(-Bip*Bj+Bi*Bjp + 2*kappa/x*Bi*Bj)
                    else if (m1 == 2 .and. m2 == 2) then
                        hq = c**2*Bip*Bjp &
                            + Bi*Bj*((Vq-c**2)**2+c**2*(-kappa*(-kappa+1)/x**2))
                    end if
                    Am(i,j) = Am(i,j) + sum(wtq*hq*jac)
                    if (m1 == m2) then
                        Bm(i,j) = Bm(i,j) + sum(wtq*Bi*Bj*jac)
                    end if
                end do
            end do
        end do
    end do
end do
! check symmetry
print *, "Checking symmetry"
do j=1,Nb
    do i=1,j-1
        if (abs(Am(i,j)-Am(j,i)) / (max(Am(i,j), Am(j,i))+tiny(1._dp)) &
                > 1e-8_dp) then
            print *, i, j, Am(i,j)-Am(j,i), Am(i,j), Am(j,i)
            call stop_error("Am not symmetric")
        end if
        if (abs(Bm(i,j)-Bm(j,i)) > 1e-12_dp) call stop_error("Bm not symmetric")
   end do
end do
end subroutine


end module


!------------------------------------------------------------------------------

program rdirac

use types, only: dp
use rdirac_assembly, only: sfem
implicit none

integer :: Ne, p, Nq, DOFs, i, kappa
real(dp), allocatable :: eigs(:)
real(dp) :: L, a, c, Z
logical :: squared

Ne = 10
p = 50
Nq = 64
L = 20
a = 1e6
Z = 92
kappa = 2
c = 137.03599907_dp
squared = .false.
call sfem(Ne, p, Nq, L, DOFs, kappa, a, c, Z, eigs, squared)
print *, "Ne:", Ne
print *, "p:", p
print *, "Nq:", Nq
print *, "DOFs:", DOFs
do i = 1, 6
    print *, i, eigs(i)
end do


end program
