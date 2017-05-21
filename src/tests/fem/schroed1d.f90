module schroed1d_assembly
use types, only: dp
use sorting, only: sort
use constants, only: pi
use feutils, only: get_parent_nodes, get_parent_quad_pts_wts, phih, dphih, &
    get_nodes, define_connect, c2fullc => c2fullc2, get_quad_pts, fe2quad
use utils, only: assert
use splines, only: iixmin, spline3pars, poly3, dpoly3
use mesh, only: meshexp
use linalg, only: eigh, eigvals
implicit none
private
public assemble_1d, assemble_1d_enr, sfem_periodic, sfem_periodic_enr, &
    sfem_non_periodic

contains

subroutine sfem_non_periodic(Ne, p, Nq, L, Nb, eigs)
integer, intent(in) :: Ne, p, Nq
real(dp), intent(in) :: L
integer, intent(out) :: Nb
real(dp), allocatable, intent(out) :: eigs(:)

integer :: Nn
! xe(i) is the 'x' coordinate of the i-th mesh node
real(dp), allocatable :: xe(:)
real(dp), allocatable :: xin(:), xiq(:), wtq(:), A(:, :), B(:, :), c(:, :), &
    phihq(:, :), dphihq(:, :), Vq(:,:), xn(:), &
    fullc(:), enrq(:,:,:), denrq(:,:,:), phipuq(:,:), dphipuq(:,:), xinpu(:), &
    uq(:,:)
integer, allocatable :: ib(:, :), in(:, :), ibenr(:,:,:)
real(dp) :: rc
integer :: i, j, iqx, u, Nenr, emin, emax

Nn = Ne*p+1
allocate(xe(Ne+1))
xe = meshexp(-10._dp, 10._dp, 1._dp, Ne) ! uniform mesh on [0, L]

allocate(xin(p+1), Vq(Nq,Ne), uq(Nq,Ne))
allocate(xinpu(2)) ! linear functions for PU
call get_parent_nodes(2, p, xin)
call get_parent_nodes(2, size(xinpu)-1, xinpu)
allocate(xiq(Nq), wtq(Nq))
call get_parent_quad_pts_wts(2, Nq, xiq, wtq)
allocate(xn(Nn))
call get_nodes(xe, xin, xn)
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

allocate(A(Nb, Nb), B(Nb, Nb), c(Nb, Nb), eigs(Nb))
allocate(fullc(Nn))

call load_potential(xe, xiq, .false., Vq)

call assemble_1d(xin, xe, ib, xiq, wtq, phihq, dphihq, Vq, A, B)
call eigh(A, B, eigs, c)
open(newunit=u, file="enrichment.txt", status="replace")
write(u, *) size(xn)
call get_quad_pts(xe,xiq,uq)
write(u, *) uq(:Nq-1,:), uq(Nq,Ne)
write(u, *) Vq(:Nq-1,:), Vq(Nq,Ne)
do i = 1, min(Nb, 20)
    call c2fullc(in, ib, c(:,i), fullc)
    if (fullc(2) < 0) fullc = -fullc
    ! Multiply by the cutoff function
    if (i == 1) then
        rc = 2.5_dp
    else if (i == 2) then
        rc = 2.5_dp
    else
        rc = 2.5_dp
    end if
    fullc = fullc*h(abs(xn-L/2), rc)
    call fe2quad(xe, xin, xiq, in, fullc, uq)
    write(u, *) uq(:Nq-1,:), uq(Nq,Ne)
end do
close(u)
end subroutine

subroutine sfem_periodic(Ne, p, Nq, L, Nb, eigs)
integer, intent(in) :: Ne, p, Nq
real(dp), intent(in) :: L
integer, intent(out) :: Nb
real(dp), allocatable, intent(out) :: eigs(:)

integer :: Nn
! xe(i) is the 'x' coordinate of the i-th mesh node
real(dp), allocatable :: xe(:)
real(dp), allocatable :: xin(:), xiq(:), wtq(:), A(:, :), B(:, :), c(:, :), &
    phihq(:, :), dphihq(:, :), Vq(:,:), xn(:), &
    fullc(:), enrq(:,:,:), denrq(:,:,:), phipuq(:,:), dphipuq(:,:), xinpu(:)
integer, allocatable :: ib(:, :), in(:, :), ibenr(:,:,:)
real(dp) :: rc
integer :: i, j, iqx, u, Nenr, emin, emax

Nn = Ne*p+1
allocate(xe(Ne+1))
xe = meshexp(0._dp, L, 1._dp, Ne) ! uniform mesh on [0, L]

allocate(xin(p+1), Vq(Nq,Ne))
allocate(xinpu(2)) ! linear functions for PU
call get_parent_nodes(2, p, xin)
call get_parent_nodes(2, size(xinpu)-1, xinpu)
allocate(xiq(Nq), wtq(Nq))
call get_parent_quad_pts_wts(2, Nq, xiq, wtq)
allocate(xn(Nn))
call get_nodes(xe, xin, xn)
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
call define_connect(3,3,Ne,p,in,ib)

Nb = maxval(ib)

allocate(A(Nb, Nb), B(Nb, Nb), c(Nb, Nb), eigs(Nb))
allocate(fullc(Nn))

call load_potential(xe, xiq, .true., Vq)

call assemble_1d(xin, xe, ib, xiq, wtq, phihq, dphihq, Vq, A, B)
eigs = eigvals(B)
print *, "cond B:", maxval(abs(eigs))/minval(abs(eigs))
call eigh(A, B, eigs, c)
end subroutine

subroutine sfem_periodic_enr(Ne, p, Nq, L, Nb, eigs)
integer, intent(in) :: Ne, p, Nq
real(dp), intent(in) :: L
integer, intent(out) :: Nb
real(dp), allocatable, intent(out) :: eigs(:)

integer :: Nn
! xe(i) is the 'x' coordinate of the i-th mesh node
real(dp), allocatable :: xe(:)
real(dp), allocatable :: xin(:), xiq(:), wtq(:), A(:, :), B(:, :), c(:, :), &
    phihq(:, :), dphihq(:, :), Vq(:,:), xn(:), &
    fullc(:), enrq(:,:,:), denrq(:,:,:), phipuq(:,:), dphipuq(:,:), xinpu(:)
integer, allocatable :: ib(:, :), in(:, :), ibenr(:,:,:)
integer :: i, j, iqx, u, Nenr, emin, emax

Nn = Ne*p+1
allocate(xe(Ne+1))
xe = meshexp(0._dp, L, 1._dp, Ne) ! uniform mesh on [0, L]

allocate(xin(p+1), Vq(Nq,Ne))
allocate(xinpu(2)) ! linear functions for PU
call get_parent_nodes(2, p, xin)
call get_parent_nodes(2, size(xinpu)-1, xinpu)
allocate(xiq(Nq), wtq(Nq))
call get_parent_quad_pts_wts(2, Nq, xiq, wtq)
allocate(xn(Nn))
call get_nodes(xe, xin, xn)
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
call define_connect(3,3,Ne,p,in,ib)
Nb = maxval(ib)

call load_potential(xe, xiq, .true., Vq)
Nenr = 2
allocate(ibenr(2,Nenr,Ne))
emin = 1
emax = Ne
!emin = 1
!emax = 2
call define_connect_enr(emin, emax, size(xinpu)-1, Nenr, Nb, ibenr)
Nb = maxval(ibenr)
allocate(enrq(Nq,Ne,Nenr))
allocate(denrq(Nq,Ne,Nenr))
call load_enrichment(xe, xiq, enrq, denrq)

allocate(A(Nb, Nb), B(Nb, Nb), c(Nb, Nb), eigs(Nb))
allocate(fullc(Nn))

call assemble_1d_enr(xin, xe, ib, ibenr, xiq, wtq, phihq, dphihq, phipuq, &
    dphipuq, Vq, enrq, denrq, A, B)
eigs = eigvals(B)
print *, "cond B:", maxval(abs(eigs))/minval(abs(eigs))
call eigh(A, B, eigs, c)
!open(newunit=u, file="wfn.txt", status="replace")
!write(u, *) xn
!do i = 1, min(Nb, 20)
!    print "(i4, f20.12)", i, eigs(i)
!    call c2fullc(in, ib, c(:,i), fullc)
!    if (fullc(2) < 0) fullc = -fullc
!    write(u, *) fullc
!end do
!close(u)
end subroutine

elemental function h(r, rc)
! C^3 cutoff function h(r)
real(dp), intent(in) :: r  ! Point at which to evaluate, must be r >= 0
real(dp), intent(in) :: rc ! Cutoff radius; h(r) = 0 for r >= rc
real(dp) :: h
if (r < rc) then
    h = 1 + 20*(r/rc)**7-70*(r/rc)**6+84*(r/rc)**5-35*(r/rc)**4
else
    h = 0
end if
end function

subroutine load_potential(xe, xiq, periodic, Vq)
real(dp), intent(in) :: xe(:), xiq(:)
logical, intent(in) :: periodic
real(dp), intent(out) :: Vq(:,:)
real(dp), allocatable :: Xn(:), Vn(:), c(:,:)
real(dp) :: jacx, x(size(xiq))
integer :: u, n, e, ip, iqx
ip = 0
do e = 1, size(xe)-1
    jacx=(xe(e+1)-xe(e))/2;
    x = xe(e) + (xiq + 1) * jacx
    do iqx = 1, size(xiq)
        Vq(iqx, e) = x(iqx)**2/2
    end do
end do
end subroutine

subroutine load_enrichment(xe, xiq, enrq, denrq)
real(dp), intent(in) :: xe(:), xiq(:)
! enrq(i,j,k) i-th quad point, j-th element, k-th enrichment
real(dp), intent(out) :: enrq(:,:,:), denrq(:,:,:)
real(dp), allocatable :: Xn(:), fn(:), c(:,:)
real(dp) :: jacx, x(size(xiq))
integer :: u, n, e, ip, i, Nenr, iqx
Nenr = size(enrq,3)
open(newunit=u, file="enrichment.txt", status="old")
read(u, *) n
allocate(Xn(n), fn(n))
allocate(c(0:4, n-1))
read(u, *) Xn
do i = 1, Nenr
    read(u, *) fn ! Load the enrichment function
    ! Interpolate using cubic splines
    call spline3pars(Xn, fn, [2, 2], [0._dp, 0._dp], c)
    ip = 0
    do e = 1, size(xe)-1
        jacx=(xe(e+1)-xe(e))/2;
        x = xe(e) + (xiq + 1) * jacx
        do iqx = 1, size(xiq)
            ip = iixmin(x(iqx), Xn, ip)
            enrq (iqx, e, i) =  poly3(x(iqx), c(:, ip))
            denrq(iqx, e, i) = dpoly3(x(iqx), c(:, ip))
        end do
    end do
end do
close(u)
end subroutine

subroutine define_connect_enr(emin, emax, p, Nenr, Nb, ibenr)
integer, intent(in) :: emin, emax, p, Nenr, Nb
integer, intent(out) :: ibenr(:,:,:)
integer :: e, i, alpha, idx
! construct nodal connectivity matrix
ibenr = 0
idx = Nb
do alpha=1,Nenr
do e=emin,emax
do i=1,p+1
    if (i == 1 .and. e == emin) cycle
    if (i == p+1 .and. e == emax) cycle
    if (i > 1) idx = idx + 1
    ibenr(i,alpha,e) = idx
end do
end do
end do
end subroutine

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


subroutine assemble_1d_enr(xin, nodes, ib, ibenr, xiq, wtq, phihq, dphihq, &
        phipuq, dphipuq, Vq, enrq, denrq, Am, Bm)
! Assemble on a 1D uniform mesh
real(dp), intent(in):: xin(:), nodes(:), xiq(:), wtq(:), &
    phihq(:, :), dphihq(:, :), Vq(:,:), enrq(:,:,:), denrq(:,:,:), &
    phipuq(:,:), dphipuq(:,:)
integer, intent(in):: ib(:, :), ibenr(:,:,:)
real(dp), intent(out):: Am(:,:), Bm(:, :)
real(dp), dimension(size(xiq), &
    size(xin)) :: phi_v, phi_dx
real(dp), dimension(size(xiq)) :: x, xp
integer :: Ne, p, e, i, j, iqx
real(dp) :: lx
integer :: ax, bx, aalpha, balpha, Nenr, u
real(dp) :: jacx, jac_det

Nenr = size(enrq,3)
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
        ! FEM x FEM
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
        ! Enrichment x FEM
        do aalpha = 1, Nenr
        do ax = 1, size(phipuq,2)
            i = ibenr(ax, aalpha, e)
            if (i == 0) cycle
            call assert(j < i)
            Am(i,j) = Am(i,j) + sum(( &
                (dphipuq(:, ax)/jacx* enrq(:,e,aalpha) &
                 +phipuq(:, ax)*denrq(:,e,aalpha)) * phi_dx(:, bx) &
                * jac_det * wtq)) / 2
            Am(i,j) = Am(i,j) + sum((Vq(:,e) * &
                phipuq(:, ax)*enrq(:,e,aalpha) * phi_v(:, bx) &
                * jac_det * wtq))
            Bm(i,j) = Bm(i,j) + sum(( &
                phipuq(:, ax)*enrq(:,e,aalpha) * phi_v(:, bx) &
                * jac_det * wtq))
        end do
        end do
    end do
    ! Enrichment x Enrichment
    do aalpha = 1, Nenr
    do ax = 1, size(phipuq,2)
        i = ibenr(ax, aalpha, e)
        if (i == 0) cycle
        do balpha = 1, Nenr
        do bx = 1, size(phipuq,2)
            j = ibenr(bx, balpha, e)
            if (j == 0) cycle
            if (j > i) cycle
            Am(i,j) = Am(i,j) + sum(( &
                (dphipuq(:, ax)/jacx* enrq(:,e,aalpha) &
                 +phipuq(:, ax)*denrq(:,e,aalpha)) &
                *(dphipuq(:, bx)/jacx* enrq(:,e,balpha) &
                 +phipuq(:, bx)*denrq(:,e,balpha)) &
                * jac_det * wtq))/2
            Am(i,j) = Am(i,j) + sum((Vq(:,e) * &
                phipuq(:, ax)*enrq(:,e,aalpha) &
                * phipuq(:, bx)*enrq(:,e,balpha) &
                * jac_det * wtq))
            Bm(i,j) = Bm(i,j) + sum(( &
                phipuq(:, ax)*enrq(:,e,aalpha) &
                * phipuq(:, bx)*enrq(:,e,balpha) &
                * jac_det * wtq))
        end do
        end do
    end do
    end do
end do
do j = 1, size(Am, 2)
    do i = 1, j-1
        Am(i, j) = Am(j, i)
        Bm(i, j) = Bm(j, i)
    end do
end do

!open(newunit=u, file="B.txt", status="replace")
!do i = 1, size(Am,1)
!    write(u,*) Bm(i, :)
!end do
!close(u)
!open(newunit=u, file="A.txt", status="replace")
!do i = 1, size(Am,1)
!    write(u,*) Am(i, :)
!end do
!close(u)
end subroutine


end module


!------------------------------------------------------------------------------

program schroed1d

use types, only: dp
use schroed1d_assembly, only: sfem_non_periodic, sfem_periodic_enr, &
    sfem_periodic
implicit none

integer :: Ne, p, Nq, DOFs, i, u, j
real(dp), allocatable :: eigs(:)
real(dp) :: L

Ne = 8
p = 50
!Nq = p+1
Nq = 64
L = 5
open(newunit=u, file="fem3.txt", status="replace")
Ne = 3
do j = 1, 8
    p = 3
    Nq = 64
    L = 20
    call sfem_non_periodic(Ne, p, Nq, L, DOFs, eigs)
    print *, "Non-periodic"
    print *, "Ne:", Ne
    print *, "p:", p
    print *, "Nq:", Nq
    print *, "DOFs:", DOFs
    do i = 1, 6
        print *, i, eigs(i)
    end do
    write(u,*) DOFs, eigs(:6)
    Ne = Ne * 2
end do
close(u)
stop "OK"

open(newunit=u, file="fem6.txt", status="replace")
do Ne = 1, 6
    p = 6
    Nq = 64
    L = 5
    call sfem_periodic(Ne, p, Nq, L, DOFs, eigs)
    print *, "Periodic"
    print *, "Ne:", Ne
    print *, "p:", p
    print *, "Nq:", Nq
    print *, "DOFs:", DOFs
    print *, "cond:", maxval(abs(eigs))/minval(abs(eigs))
    do i = 1, 2
        print *, i, eigs(i)
    end do
    write(u,*) DOFs, p, Ne, Nq, L, eigs(:2)
end do
close(u)

open(newunit=u, file="fem3.txt", status="replace")
do Ne = 2, 12, 2
    p = 3
    Nq = 64
    L = 5
    call sfem_periodic(Ne, p, Nq, L, DOFs, eigs)
    print *, "Periodic"
    print *, "Ne:", Ne
    print *, "p:", p
    print *, "Nq:", Nq
    print *, "DOFs:", DOFs
    print *, "cond:", maxval(abs(eigs))/minval(abs(eigs))
    do i = 1, 2
        print *, i, eigs(i)
    end do
    write(u,*) DOFs, p, Ne, Nq, L, eigs(:2)
end do
close(u)


open(newunit=u, file="sfem2.txt", status="replace")
do p = 1, 63
    Ne = 2
    Nq = 64
    L = 5
    call sfem_periodic(Ne, p, Nq, L, DOFs, eigs)
    print *, "Periodic"
    print *, "Ne:", Ne
    print *, "p:", p
    print *, "Nq:", Nq
    print *, "DOFs:", DOFs
    print *, "cond:", maxval(abs(eigs))/minval(abs(eigs))
    do i = 1, 2
        print *, i, eigs(i)
    end do
    write(u,*) DOFs, p, Ne, Nq, L, eigs(:2)
end do
close(u)

open(newunit=u, file="sfem_enr.txt", status="replace")
do p = 1, 15
    Ne = 2
    Nq = 64
    L = 5
    call sfem_periodic_enr(Ne, p, Nq, L, DOFs, eigs)
    print *, "Periodic + enrichment"
    print *, "Ne:", Ne
    print *, "p:", p
    print *, "Nq:", Nq
    print *, "DOFs:", DOFs
    print *, "cond:", maxval(abs(eigs))/minval(abs(eigs))
    do i = 1, 2
        print *, i, eigs(i)
    end do
    write(u,*) DOFs, p, Ne, Nq, L, eigs(:2)
end do

end program
