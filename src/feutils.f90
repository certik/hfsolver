module feutils
! Various Finite Element (FE) utilities
! Based on code written by John E. Pask, LLNL.
use types
use utils, only: stop_error
use quadrature, only: gauss_pts, gauss_wts, lobatto_wts, lobatto_pts
use solvers, only: solve_sym
implicit none
private
public get_nodes, define_connect, get_quad_pts, define_connect_n, &
    define_connect_np, get_parent_quad_pts_wts, &
    fe2quad, c2fullc, c2fullc2, fe_evaluate, calc_proj_error, &
    get_parent_nodes, phih, dphih

contains

subroutine get_parent_quad_pts_wts(qtype, Nq, xiq, wtq)
! Quadrature points and weights
integer, intent(in) :: qtype      ! type of quadrature: 1= Gauss, 2= Lobatto
integer, intent(in) :: Nq         ! number or quadrature points/weights
real(dp), intent(out) :: xiq(:)   ! quadrature points
real(dp), intent(out) :: wtq(:)   ! quadrature weights
select case(qtype)
    case(1)  ! Gauss
        xiq = gauss_pts(Nq)
        wtq = gauss_wts(Nq)
    case(2)  ! Lobatto
        xiq = lobatto_pts(Nq)
        wtq = lobatto_wts(Nq)
    case default
        call stop_error("Invalid quadrature type.")
end select
end subroutine

subroutine get_nodes(xe,xin,xn)
! generates basis nodes xn via affine mapping of parent nodes xin to elements xe
real(dp), intent(in) :: xe(:)     ! elements: xe(i/i+1) = coord of left/right boundary of ith element
real(dp), intent(in) :: xin(:)    ! parent basis nodes: xin(i) = coord of ith parent basis fn node
real(dp), intent(out) :: xn(:)    ! basis nodes: xn(i) = coordinate of ith basis fn node
real(dp) :: xa,xb                ! left and right element boundaries
integer Ne,p,i,j
Ne=size(xe)-1
p=size(xin)-1
xn(1)=xe(1) ! left-most node
do i=1,Ne   ! for each element ...
    xa=xe(i); xb=xe(i+1)
    xn(i*p+1)=xb   ! right-end node
    do j=2,p
        xn((i-1)*p+j)=(xin(j)+1)/2*(xb-xa)+xa ! internal nodes: affine mapping
    end do
end do
end subroutine

!-------------------------------------------------------------------------------------------------!

subroutine define_connect(bca,bcb,Ne,p,in,ib)
! constructs connectivity matrices in and ib defining local-global node and basis
! correspondence, respectively
integer, intent(in):: bca        ! boundary condition at x=a: 1=Dirichlet, 2=Neumann,
   ! 3=periodic, 4=antiperiodic. If bca=3 or 4, bcb must be set correspondingly
integer, intent(in):: bcb        ! boundary condition at x=b: 1=Dirichlet, 2=Neumann,
   ! 3=periodic, 4=antiperiodic. If bcb = 3 or 4, bca must be set correspondingly
integer, intent(in):: Ne         ! number of elements
integer, intent(in):: p          ! order of FE/SE basis
integer, intent(out):: in(:,:)   ! nodal connectivity: in(i,j) = index of basis node
                                 ! corresponding to local node i of element j
integer, intent(out):: ib(:,:)  ! basis connectivity: ib(i,j) = index of basis function
   ! associated with local basis function i of element j. 0 = no associated basis fn.
   ! -1 = associated with antiperiodic basis fn, with negative multiplier
integer i,e
! check boundary condition consistency
if ((bca>2 .or. bcb>2) .and. bcb/=bca) then
   write(*,'(1x,a,i0,a,i0)') "Error: bca = ", bca, " /= bcb = ", bcb
   call stop_error("stop")
end if
! construct nodal connectivity matrix
do e=1,Ne
   do i=1,p+1
      in(i,e)=(e-1)*p+i
   end do
end do
! construct basis connectivity matrix
   ! construct for natural BCs
ib=in
   ! modify for other BCs
if (bca==1) ib=ib-1        ! Dirichlet at x=a -> omit basis fn at node 1
select case(bcb)           ! x=b....
   case(1); ib(p+1,Ne)=0     ! Dirichlet
   case(3); ib(p+1,Ne)=1     ! periodic
   case(4); ib(p+1,Ne)=-1    ! antiperiodic
end select
end subroutine

subroutine define_connect_n(bca, bcb, Ne, p, Nm, in, ib)
! constructs connectivity matrices in and ib defining local-global node and
! basis correspondence, respectively. Works for arbitrary number of meshes
integer, intent(in):: bca(:)        ! bca(m) boundary condition at x=a for mesh
    ! m: 1=Dirichlet, 2=Neumann
integer, intent(in):: bcb(:)        ! bca(m) boundary condition at x=b for mesh
    ! m: 1=Dirichlet, 2=Neumann
integer, intent(in):: Ne         ! number of elements
integer, intent(in):: p          ! order of FE/SE basis
integer, intent(in):: Nm         ! number of meshes
integer, intent(out):: in(:,:,:)   ! nodal connectivity: in(i,j,k) = index of
    ! basis node corresponding to local node i of element j of mesh k
integer, intent(out):: ib(:,:,:)  ! basis connectivity: ib(i,j,k) = index of
    ! basis function associated with local basis function i of element j of
    ! mesh k.  0 = no associated basis fn
integer :: i, e, m, dofs
! construct nodal connectivity matrix
dofs = 0
do m = 1, Nm
    do e = 1, Ne
        do i = 1, p+1
            dofs = dofs + 1
            in(i, e, m) = dofs
        end do
        dofs = dofs - 1 ! connect the shape functions across the element
    end do
    dofs = dofs + 1 ! don't connect the shape functions across meshes
end do
! construct basis connectivity matrix
dofs = 0
do m = 1, Nm
    do e = 1, Ne
        do i = 1, p+1
            if (i == 1   .and. e == 1  .and. bca(m) == 1) then
                ib(i, e, m) = 0 ! Dirichlet
                cycle
            end if
            if (i == p+1 .and. e == Ne .and. bcb(m) == 1) then
                ib(i, e, m) = 0 ! Dirichlet
                cycle
            end if
            dofs = dofs + 1
            ib(i, e, m) = dofs
        end do
        dofs = dofs - 1 ! connect the shape functions across the element
    end do
    dofs = dofs + 1 ! don't connect the shape functions across meshes
end do
end subroutine

subroutine define_connect_np(bca, bcb, Ne, p, Nm, in, ib)
! constructs connectivity matrices in and ib defining local-global node and
! basis correspondence, respectively. Works for arbitrary number of meshes
integer, intent(in):: bca(:)        ! bca(m) boundary condition at x=a for mesh
    ! m: 1=Dirichlet, 2=Neumann
integer, intent(in):: bcb(:)        ! bca(m) boundary condition at x=b for mesh
    ! m: 1=Dirichlet, 2=Neumann
integer, intent(in):: Ne         ! number of elements
integer, intent(in):: p(:)          ! order of FE/SE basis
integer, intent(in):: Nm         ! number of meshes
integer, intent(out):: in(:,:,:)   ! nodal connectivity: in(i,j,k) = index of
    ! basis node corresponding to local node i of element j of mesh k
integer, intent(out):: ib(:,:,:)  ! basis connectivity: ib(i,j,k) = index of
    ! basis function associated with local basis function i of element j of
    ! mesh k.  0 = no associated basis fn
integer :: i, e, m, dofs
! Initialize the arrays to -1
in = -1
ib = -1
! construct nodal connectivity matrix
dofs = 0
do m = 1, Nm
    do e = 1, Ne
        do i = 1, p(m)+1
            dofs = dofs + 1
            in(i, e, m) = dofs
        end do
        dofs = dofs - 1 ! connect the shape functions across the element
    end do
    dofs = dofs + 1 ! don't connect the shape functions across meshes
end do
! construct basis connectivity matrix
dofs = 0
do m = 1, Nm
    do e = 1, Ne
        do i = 1, p(m)+1
            if (i == 1   .and. e == 1  .and. bca(m) == 1) then
                ib(i, e, m) = 0 ! Dirichlet
                cycle
            end if
            if (i == p(m)+1 .and. e == Ne .and. bcb(m) == 1) then
                ib(i, e, m) = 0 ! Dirichlet
                cycle
            end if
            dofs = dofs + 1
            ib(i, e, m) = dofs
        end do
        dofs = dofs - 1 ! connect the shape functions across the element
    end do
    dofs = dofs + 1 ! don't connect the shape functions across meshes
end do
end subroutine

subroutine get_quad_pts(xe,xiq,xq)
! generates quadrature points xq via affine mapping of parent quad points xiq to elements xe
real(dp), intent(in):: xe(:)     ! elements: xe(i/i+1) = coord of left/right boundary of ith element
real(dp), intent(in):: xiq(:)    ! parent quad pts
real(dp), intent(out):: xq(:,:)  ! quad pts: xq(i,j) = coordinate of ith point in jth element
real(dp) xa,xb                ! left and right element boundaries
integer ie
do ie=1,size(xe)-1
   xa=xe(ie); xb=xe(ie+1)
   xq(:,ie)=(xb-xa)/2*xiq+(xb+xa)/2 ! affine transformation
end do
end subroutine

subroutine fe2quad(xe, xin, xiq, in, fullu, uq)
! transforms fullu from FE-coefficient to quadrature-grid representation.
! fullu is a full FE coefficient vector, having values for all nodes in the mesh,
! including domain-boundary nodes.
real(dp), intent(in) :: xe(:)     ! elements: xe(i/i+1) = coord of left/right boundary of ith element
real(dp), intent(in) :: xin(:)    ! parent basis nodes: xin(i) = coordinate of ith parent basis node
real(dp), intent(in) :: xiq(:)    ! quadrature points
integer, intent(in) :: in(:,:)    ! nodal connectivity: in(i,j) = index of basis node
                                 ! corresponding to local node i of element j
real(dp), intent(in) :: fullu(:)  ! FE coefficient representatin of fullu: full vector, including
   ! values for domain-boundary nodes / basis functions
real(dp), intent(out) :: uq(:,:)  ! quadrature-grid representatin of fullu
   ! uq(i,j) = value at ith quadrature point of jth element
real(dp) :: phihq(size(xiq), size(xin)) ! parent basis fn values at quadrature points:
   ! phihq(i,j) = value of jth function at ith quadrature point
integer :: ie, iln, iq                ! element, local node, quad point indices

! tabulate parent basis at quadrature points
do iln = 1, size(xin)
    do iq = 1, size(xiq)
        phihq(iq, iln) = phih(xin, iln, xiq(iq))
    end do
end do
! evaluate at quad points in each element
do ie = 1, size(xe)-1
    uq(:, ie) = 0
    do iln = 1, size(xin)
        uq(:, ie) = uq(:, ie) + fullu(in(iln, ie)) * phihq(:, iln)
    end do
end do
end subroutine

subroutine c2fullc(in, ib, c, fullc)
! Converts FE coefficient vector to full coefficient vector
integer, intent(in) :: in(:,:,:)   ! nodal connectivity: in(i,j) = index of basis
   ! node
                                 ! corresponding to local node i of element j
integer, intent(in) :: ib(:,:,:)   ! basis connectivity: ib(i,j) = index of basis
   ! function associated with local basis function i of element j. 0 = no
   ! associated basis fn.
real(dp), intent(in) :: c(:) ! coefficient vector with regards to ib
real(dp), intent(out) :: fullc(:) ! full coefficients vector with regards to in
integer :: m, e, i
do m = 1, size(in, 3)
    do e = 1, size(in, 2)
        do i = 1, size(in, 1)
            if (in(i, e, m) == -1) cycle
            if (ib(i, e, m) == 0) then
                fullc(in(i, e, m)) = 0 ! Dirichlet
            else
                fullc(in(i, e, m)) = c(ib(i, e, m))
            end if
        end do
    end do
end do
end subroutine

subroutine c2fullc2(in, ib, c, fullc)
! Converts FE coefficient vector to full coefficient vector
! It puts 0 for Dirichlet boundary conditions (ib==0), otherwise it just copies
! the coefficients.
integer, intent(in) :: in(:,:)   ! nodal connectivity: in(i,j) = index of basis
   ! node
                                 ! corresponding to local node i of element j
integer, intent(in) :: ib(:,:)   ! basis connectivity: ib(i,j) = index of basis
   ! function associated with local basis function i of element j. 0 = no
   ! associated basis fn.
real(dp), intent(in) :: c(:) ! coefficient vector with regards to ib
real(dp), intent(out) :: fullc(:) ! full coefficients vector with regards to in
integer :: e, i
do e = 1, size(in, 2)
    do i = 1, size(in, 1)
        if (ib(i, e) == 0) then
            fullc(in(i, e)) = 0 ! Dirichlet
        else
            fullc(in(i, e)) = c(ib(i, e))
        end if
    end do
end do
end subroutine

subroutine fe_evaluate(xe, xin, xn, in, fullu, xout, yout)
! Evaluates the FE solution (fullu) on a grid.
! fullu is a full FE coefficient vector, having values for all nodes in the mesh,
! including domain-boundary nodes.
real(dp), intent(in):: xe(:)     ! elements: xe(i/i+1) = coord of left/right boundary of ith element
real(dp), intent(in):: xin(:)    ! parent basis nodes: xin(i) = coordinate of ith parent basis node
real(dp), intent(in):: xn(:)     ! basis nodes: xn(i) = coordinate of ith basis node
integer, intent(in):: in(:,:)    ! nodal connectivity: in(i,j) = index of basis node
                                 ! corresponding to local node i of element j
real(dp), intent(in):: fullu(:, :)  ! FE coefficient representatin of fullu: full vector, including
   ! values for domain-boundary nodes / basis functions
real(dp), intent(in):: xout(:)   ! grid points to evaluate the solution at
real(dp), intent(out):: yout(:, :)  ! values of the solution at 'xout'
real(dp) :: phihx(size(xin)) ! parent basis fn values at a point
integer :: p, n, e, i, j
real(dp) :: x, xa, xb, xi

p = size(xin) - 1
e = 1
do i = 1, size(xout)
    ! get x and element containing x
    x = xout(i)
    call getElement(xn, p, x, e)
    ! get parent coordinates xi
    xa = xe(e)
    xb = xe(e+1)
    xi = (x-xa)/(xb-xa)*2-1
    ! get parent basis values at xi
    do n= 1, p+1
        phihx(n) = phih(xin, n, xi)
    end do
    ! get eigenfunction value at x
    do j = 1, size(fullu, 2)
        yout(i, j) = sum(fullu(in(:,e), j) * phihx)
    end do
end do
end subroutine

subroutine getElement(xn,p,x,e)
! returns element containing point x, at or after element e
real(dp), intent(in):: xn(:)  ! basis nodes: xn(i) = coordinate of ith basis node
integer, intent(in):: p       ! order of FE/SE basis
real(dp), intent(in):: x      ! point x
integer, intent(inout):: e    ! input: index of element to start search
                              ! output: index of element containing x, at or after element e
integer Ne                 ! number of elements
Ne=(size(xn)-1)/p
if (e<1 .or. e>Ne) e=1
do
   if (x>=xn((e-1)*p+1) .and. x<=xn(e*p+1)) exit
   e=e+1
   if (e>Ne) call stop_error("getElement error: x not in mesh.")
end do
end subroutine

function calc_proj_error(p, xe, xiq, wtq, uq) result(error)
! Calculates the L2 error of the projected solution on the mesh "p"
integer, intent(in) :: p
real(dp), intent(in):: xe(:)     ! elements: xe(i/i+1) = coord of left/right boundary of ith element
real(dp), intent(in):: xiq(:), wtq(:)    ! quadrature points and weights
real(dp), intent(in):: uq(:,:)  ! FE coefficient representatin of fullu: full vector, including
real(dp) :: error(size(xe)-1) ! errors for each element
   ! values for domain-boundary nodes / basis functions
real(dp) phihq(size(xiq),p+1) ! parent basis fn values at quadrature points: 
   ! phihq(i,j) = value of jth function at ith quadrature point
integer :: ie,iln,iq                 ! element, local node, quad point indices
integer :: i, j, Ne, dofs
real(dp) :: xin(p+1), al, be, jac, xa, xb, norm
real(dp), allocatable :: A(:, :), f(:), x(:), xn(:)
integer, allocatable :: in(:, :), ib(:, :)
real(dp) :: projq(size(uq, 1), size(uq, 2))

if (.not. p >= 1) call stop_error("p >= 1 required")

! define parent basis
call get_parent_nodes(2, p, xin)

Ne = size(xe) - 1

! define basis
allocate(xn(Ne*p+1))
call get_nodes(xe, xin, xn)

allocate(in(p+1,Ne),ib(p+1,Ne))
call define_connect(2, 2, Ne, p, in, ib)

! tabulate parent basis at quadrature points
do iln = 1, p+1
   do iq = 1, size(xiq)
      phihq(iq,iln) = phih(xin,iln,xiq(iq))
   end do
end do

dofs = maxval(in)
!print *, "DOFS=", dofs
allocate(A(dofs, dofs), f(dofs), x(dofs))

A = 0
f = 0
do ie = 1, Ne ! sum over elements
    xa=xe(ie); xb=xe(ie+1)
    al=(xb-xa)/2; be=(xb+xa)/2
    jac=al
    do i = 1, p+1
        do j = 1, p+1
            A(in(i,ie), in(j,ie)) = A(in(i,ie), in(j,ie)) + &
                sum(wtq*phihq(:, i)*phihq(:, j)*jac)
        end do
        f(in(i,ie)) = f(in(i,ie)) + sum(wtq*phihq(:, i)*uq(:, ie)*jac)
    end do
end do
! Solve A x = f
x = solve_sym(A, f)
call fe2quad(xe,xin,xiq,in,x,projq)
norm = 0
do ie = 1, Ne ! sum over elements
    xa=xe(ie); xb=xe(ie+1)
    al=(xb-xa)/2; be=(xb+xa)/2
    jac=al
    !print *, "el =", ie, sqrt(sum(wtq * abs(uq(:, ie)-projq(:, ie))**2 * jac))
    error(ie) = sum(wtq * abs(uq(:, ie)-projq(:, ie))**2 * jac)
    norm = norm + sum(wtq * abs(uq(:, ie))**2 * jac)
end do
norm = sqrt(norm)
error = sqrt(error)/norm
!open(newunit=ff, file="proj.dat", status="replace")
!write(ff, *) xn
!write(ff, *) x
!close(ff)
end function

subroutine get_parent_nodes(btype, p, xin)
integer, intent(in):: btype    ! type of basis: 1=uniform node, 2=Lobatto node
integer, intent(in):: p        ! order of FE/SE basis
! parent basis nodes: xin(i) = coordinate of ith parent basis node:
real(dp), intent(out):: xin(:)
integer i
if (p < 1) call stop_error("Error: p < 1.")
if (size(xin) /= p+1) call stop_error("Error: size(xin) /= p+1")
select case (btype)
    case(1)  ! uniform
        do i = 1, p+1
            xin(i) = real(i-1, dp)/p*2-1
        end do
    case(2)  ! Lobatto
        xin=lobatto_pts(p+1)
    case default
        call stop_error("Error: invalid basis type.")
end select
end subroutine

real(dp) function phih(xin, n, xi)
! "phi hat": nth Lagrange polynomial with nodes xin at point xi in [-1,1]
real(dp), intent(in) :: xin(:) ! polynomial nodes
integer, intent(in) :: n       ! polynomial index
real(dp), intent(in) :: xi     ! point in [-1,1] at which to evaluate phih
integer :: i
! compute nth polynomal: 1 at node n, 0 at all others
phih = 1
do i = 1, size(xin)
    if (i == n) cycle
    phih = phih * (xi-xin(i))/(xin(n)-xin(i))
end do
end function

real(dp) function dphih(xin, n, xi)
! "d phi hat": derivative of nth Lagrange polynomial with nodes xin at point xi
! in [-1,1]
real(dp), intent(in) :: xin(:) ! polynomial nodes
integer, intent(in) :: n       ! polynomial index
real(dp), intent(in) :: xi     ! point in [-1,1] at which to evaluate dphih
real(dp) :: term
real(dp) :: tmp(size(xin))
integer :: i, j
do i = 1, size(xin)
    if (i==n) cycle
    tmp(i) = (xi-xin(i))/(xin(n)-xin(i))
end do
! compute derivative of nth polynomal
dphih = 0
do j = 1, size(xin)
    if (j == n) cycle
    term = 1 / (xin(n)-xin(j))
    do i = 1, size(xin)
        if (i==n .or. i==j) cycle
        term = term * tmp(i)
    end do
    dphih = dphih + term
end do
end function

end module
