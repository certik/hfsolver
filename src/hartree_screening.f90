module hartree_screening
! Hartree screening function
use types, only: dp
use utils, only: stop_error
use feutils, only: define_connect, c2fullc2, fe2quad, get_parent_nodes, &
        phih, dphih
use solvers, only: solve_sym
implicit none
private
public hartree_y

contains

function hartree_y(k, p, xe, xiq, wtq, f) result(Yq)
! Calculates the Hartree screening function:
!   Y^k(f(r), r) = r*\int_0^oo dr' f(r') * min(r, r')^k / max(r, r')^(k+1)
! This function is calculated by solving an equivalent differential equation:
!   -Y^k'' + (k(k+1)/r^2) * Y^k = (2k+1)/r * f(r)
! with boundary conditions at 0 and oo (imposed by the above integral):
!   Y^k(0) = 0
!   Y^k'(oo) + (k/r) * Y^k(oo) = 0
! using FEM on the mesh 'xe', polynomial order 'p'. The function 'f' is given
! on the quadrature grid xiq/wtq. The Y^k(r) function is returned on the same
! quadrature grid. 'k' is the parameter of Y^k.
!
! Hint: use such quadrature xiq/wtq and such 'p', so that you can represent the
! function 'f' exactly. If for example f(r) = P^2(r) and P(r) is of order 'p',
! then use '2*p' for hartree_y().
!
! Note: If you want to solve the following radial Poisson equation:
!   V_H'' + (2/r)*V_H = -4*pi*n
!   (V_H*r)''         = -4*pi*n*r
! then V_H = Y^0(4*pi*n*r^2, r) / r
integer, intent(in) :: k, p
real(dp), intent(in) :: xe(:), xiq(:), wtq(:)
real(dp), intent(in) :: f(:,:)
real(dp) :: Yq(size(xiq), size(xe)-1)
real(dp), allocatable :: Am(:, :), bv(:), u(:), fullc(:), xin(:)
integer, allocatable :: ib(:, :), in(:, :)
integer :: Ne, Nb
Ne = size(xe)-1
allocate(xin(p+1))
call get_parent_nodes(2, p, xin)
allocate(in(p+1, Ne), ib(p+1, Ne))
call define_connect(1, 2, Ne, p, in, ib)
Nb = maxval(ib)
allocate(Am(Nb, Nb), bv(Nb), u(Nb))
call assemble_hartree_y(k, f, xin, xe, ib, xiq, wtq, Am, bv)
! solve
u = solve_sym(Am, bv)
! transform solution to quadrature grid
allocate(fullc(maxval(in)))
call c2fullc2(in, ib, u, fullc)
call fe2quad(xe, xin, xiq, in, fullc, Yq)
end function

subroutine assemble_hartree_y(k, f_vals, xin, xe, ib, xiq, wtq, Am, bv)
! forms system equation matrices corresponding to the problem
! -u''(r) + k(k+1)/r^2 * u(r) = (2k+1)/r * f(x)
! subject to boundary conditions consistent with basis specified by ib
! Note: u(r) is the hartree Y^k(r) screening function
! The weak formulation is:
! \int u'(r)*v'(r) + k(k+1)/r^2 * u(r) v(r) \d r =
!    =\int (2k+1)/r * f(x) * v(x) \d r
! It solves the following weak formulation
integer, intent(in) :: k            ! k parameter in the equation
real(dp), intent(in):: f_vals(:,:)   ! f(x) at quadrature points:
   ! f_vals(i,j) = value at ith point in jth element
real(dp), intent(in) :: xin(:)       ! parent basis nodes
real(dp), intent(in) :: xe(:)        ! element coordinates
integer, intent(in) :: ib(:,:)       ! basis connectivity: ib(i,j) = index of
    ! basis function associated with local basis function i of element j.
    ! 0 = no associated basis fn.
real(dp), intent(in) :: xiq(:)       ! quadrature points
real(dp), intent(in) :: wtq(:)       ! quadrature weights
real(dp), intent(out) :: Am(:,:)     ! system matrix: Am c = bv
real(dp), intent(out) :: bv(:)       ! source vector: Am c = bv
integer :: Ne, Nb                ! number of elements, basis functions
integer :: p                     ! order of FE/SE basis
integer :: e                     ! element index
integer :: i,j                   ! basis fn indices
integer :: al,be                 ! "alpha", "beta": local basis fn indices
integer :: iq                    ! quadrature point index
real(dp) :: xa,xb                ! element boundary node coordinates
real(dp) :: jac                  ! Jacobian of transformation from parent
    ! coords xi in [-1,1] to coords x in [xa,xb]: x = (xb-xa)/2 xi + (xb+xa)/2
real(dp) :: phihq(size(xiq),size(xin)), dphihq(size(xiq),size(xin))   ! parent
    ! basis fns and derivs at quadrature points:
    ! phihq(i,j) = value of jth function at ith quadrature point
real(dp), dimension(size(xiq)) :: intq, aq, bq, xq

! initializations
Ne = size(ib, 2)
Nb = maxval(ib)
p = size(xin)-1
if (size(xin) /= size(ib,1)) &
    call stop_error("Error: inconsistent parent node and connectivity dimensions.")
if (size(Am,1) /= Nb .or. size(bv,1) /= Nb) &
    call stop_error("Error: size of Am and/or bv inconsistent with Nb.")
! tabulate parent basis and derivatives at quadrature points
do al = 1, p+1
   do iq = 1, size(xiq)
      phihq(iq, al) = phih(xin, al, xiq(iq))
      dphihq(iq, al) = dphih(xin, al, xiq(iq))
   end do
end do
! accumulate Am matrix and bv vector
! compute lower triangle
Am = 0; bv = 0
do e = 1, Ne
   xa = xe(e)
   xb = xe(e+1)
   jac = (xb - xa)/2  ! affine mapping
   ! tabulate a(x), b(x) at quadrature points
   xq = (xiq(:)+1)/2*(xb-xa)+xa ! affine mapping
   aq = k*(k+1)/xq**2
   bq = (2*k+1)/xq
   ! compute matrix/vector elements (integrals transformed to [-1,1])
   do al = 1, p+1
      i = ib(al, e)
      if (i == 0) cycle              ! omit boundary basis fns for Dirichlet BCs
      bv(i) = bv(i) + sum(wtq*bq*f_vals(:, e)*phihq(:,al)*jac)
      do be = 1, p+1
         j = ib(be, e)
         if (j == 0) cycle           ! omit boundary basis fns for Dirichlet BCs
         if (j > i) cycle            ! compute only lower triangles
         intq = dphihq(:, al) * dphihq(:, be) / jac**2   &
                    + aq * phihq(:, al) * phihq(:, be)
         Am(i, j) = Am(i, j) + sum(wtq*intq*jac)
      end do
   end do
end do
! We add the surface integral term, which is equal to:
! k/r_max * phi(r_max)*phi(r_max)
! but phi(r_max)=1, so we only get k/r_max contribution to the last diagonal
! matrix element. Note: Using numerical tests, this correction is very
! important for smaller r_max (~70 atomic units), where it makes a difference
! on the 6th significant digit. For large r_max (~500 atomic units) it doesn't
! seem to mater much. For rmax->oo or k==0 this correction is 0.
i = ib(p+1, Ne)
if (i == 0) call stop_error("Inconsistent BC")
Am(i, i) = Am(i, i) + k / xe(Ne+1)
! fill in upper triangle
do j = 1, Nb
   do i = 1, j-1
      Am(i, j) = Am(j, i)
   end do
end do
end subroutine

end module
