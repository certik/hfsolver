module fe
use types, only: dp
use constants, only: pi
use utils, only: stop_error
use scf, only: ijkl2intindex
use feutils, only: get_quad_pts, phih, dphih
use utils, only: str
use hartree_screening, only: hartree_y
use special_functions, only: wigner3j
implicit none
private
public get_basis, feints

contains

subroutine get_basis(Nmax, Lmax, nlist, llist, mlist)
integer, intent(in) :: Nmax, Lmax
integer, allocatable, intent(out) :: nlist(:), llist(:), mlist(:)
integer :: i, n, l, m, a
a = Nmax * (Lmax+1)**2
allocate(nlist(a), llist(a), mlist(a))
i = 1
do n = 1, Nmax
    do l = 0, Lmax
        do m = -l, l
            nlist(i) = n
            llist(i) = l
            mlist(i) = m
            i = i + 1
        end do
    end do
end do
end subroutine

subroutine feints(Z, xin, xe, ib, xiq, wtq, S, T, V, slater)
! Calculates S, T, V arrays and slater.
! slater only has R^k(n1, n2, n3, n4) integrals, the same for all "l". As such
! the n1 index only runs over "n", not "l".
integer, intent(in) :: Z
real(dp), intent(in) :: xin(:)       ! parent basis nodes
real(dp), intent(in) :: xe(:)        ! element coordinates
integer, intent(in) :: ib(:, :)       ! basis connectivity: ib(i,j) = index of
   ! basis function associated with local basis function i of element j. 0 = no
   ! associated basis fn.
real(dp), intent(in) :: xiq(:)       ! quadrature points
real(dp), intent(in) :: wtq(:)       ! quadrature weights
real(dp), intent(out) :: S(:, :, 0:), T(:, :, 0:), V(:, :, 0:)
real(dp), intent(out) :: slater(:, 0:)
real(dp), allocatable :: slater_(:, :)
integer :: k, l, Nmax, Lmax
Nmax = maxval(ib)
Lmax = ubound(S, 3)
k = Nmax*(Nmax+1)/2
allocate(slater_(k*(k+1)/2, 0:2*Lmax))
call assemble_radial_STV(Z, Lmax, xin, xe, ib, xiq, wtq, S(:, :, 0), &
    T, V(:, :, 0))
do l = 1, Lmax
    S(:, :, l) = S(:, :, 0)
    V(:, :, l) = V(:, :, 0)
end do
call assemble_radial_int2(xin, xe, ib, xiq, wtq, ubound(slater_, 2), slater_)
slater = slater_
end subroutine

subroutine assemble_radial_STV(Z, Lmax, xin, xe, ib, xiq, wtq, S, T, V)
integer, intent(in) :: Z, Lmax
real(dp), intent(in) :: xin(:)       ! parent basis nodes
real(dp), intent(in) :: xe(:)        ! element coordinates
integer, intent(in) :: ib(:, :)       ! basis connectivity: ib(i,j) = index of
   ! basis function associated with local basis function i of element j. 0 = no
   ! associated basis fn.
real(dp), intent(in) :: xiq(:)       ! quadrature points
real(dp), intent(in) :: wtq(:)       ! quadrature weights
real(dp), intent(out) :: S(:, :), T(:, :, 0:), V(:, :)
integer Ne, Nb, p, e, i, j, al, be, iq, l
real(dp) xa, xb, jac
real(dp), dimension(size(xiq), size(xin)) :: phihq, dphihq
real(dp), dimension(size(xiq)) :: hq, x

Ne = size(xe)-1
Nb = maxval(ib)
p = size(xin)-1
! tabulate parent basis and derivatives at quadrature points
do al = 1, p+1
    do iq = 1, size(xiq)
        phihq(iq, al)  =  phih(xin, al, xiq(iq))
        dphihq(iq, al) = dphih(xin, al, xiq(iq))
    end do
end do
S = 0
T = 0
V = 0
do e = 1, Ne
   xa = xe(e)
   xb = xe(e+1)
   jac = (xb-xa)/2  ! affine mapping
   x = (xiq+1)/2*(xb-xa)+xa
   do be = 1, p+1
      j = ib(be, e)
      if (j==0) cycle               ! omit boundary basis fns for Dirichlet BCs
      do al = 1, p+1
         i = ib(al, e)
         if (i==0) cycle            ! omit boundary basis fns for Dirichlet BCs
         if (j>i) cycle             ! compute only lower triangles
         hq = phihq(:, al) * phihq(:,be)
         S(i, j) = S(i, j) + sum(wtq * hq * jac)
         do l = 0, Lmax
             hq = 0.5_dp * dphihq(:, al)*dphihq(:, be)/jac**2 &
                    + phihq(:, al) * l*(l+1)/(2*x**2) * phihq(:,be)
             T(i, j, l) = T(i, j, l) + sum(wtq * hq * jac)
         end do
         hq = phihq(:, al) * (-Z/x) * phihq(:, be)
         V(i, j) = V(i, j) + sum(wtq * hq * jac)
      end do
   end do
end do
! fill in upper triangles
do j = 1, Nb
    do i = 1, j-1
        S(i, j) = S(j, i)
        do l = 0, Lmax
            T(i, j, l) = T(j, i, l)
        end do
        V(i, j) = V(j, i)
    end do
end do
end subroutine

pure integer function sym_idx(i, j) result(k)
! Converts a symmetric matrix index (i, j) into a vector index (k).
! The upper triangular part of the symmetrix matrix is stored in a
! one-dimensional array, packed by columns.
! For a symmetrix n x n matrix, the vector will have n*(n+1)/2 elements.
! As such, for 1 <= i, j <= n, it will return 1 <= k <= n*(n+1)/2.
integer, intent(in) :: i, j
if (j >= i) then
    k = i + j * (j-1) / 2
else
    k = j + i * (i-1) / 2
end if
end function

real(dp) function factor(l1, l1p, l2, l2p, m1, m1p, m2, m2p, k) result(r)
integer, intent(in) :: l1, l1p, l2, l2p, m1, m1p, m2, m2p, k
integer :: q
r = 0
do q = -k, k
    r = r + (-1)**(m1+m2+q) &
        * wigner3j(l1, l1p, k, -m1, m1p, +q) &
        * wigner3j(l2, l2p, k, -m2, m2p, -q)
end do
if (abs(r-factor2(l1, l1p, l2, l2p, m1, m1p, m2, m2p, k)) > 1e-9_dp) then
    stop "FAIL"
end if
end function

real(dp) function factor2(l1, l1p, l2, l2p, m1, m1p, m2, m2p, k) result(r)
integer, intent(in) :: l1, l1p, l2, l2p, m1, m1p, m2, m2p, k
if (k >= abs(m1-m1p) .and. m1+m2-m1p-m2p == 0) then
    r = (-1)**(m1+m2p) &
        * wigner3j(l1, l1p, k, -m1, m1p, m1-m1p) &
        * wigner3j(l2, l2p, k, -m2, m2p, m2-m2p)
else
    r = 0
end if
end function

subroutine assemble_radial_int2(xin, xe, ib, xiq, wtq, max_k, int2)
! Calculates the 2 particle integral:
! \int dx dy P1(x)*P2(x)*P1(y)*P2(y)*min(x,y)^k / max(x, y)^(k+1)
real(dp), intent(in) :: xin(:), xe(:), xiq(:), wtq(:)
integer, intent(in) :: ib(:,:)       ! basis connectivity: ib(i,j) = index of
    ! basis function associated with local basis function i of element j.
    ! 0 = no associated basis fn.
integer, intent(in) :: max_k
real(dp), intent(out) :: int2(:, 0:)
integer :: p
integer :: Ne, Nb                ! number of elements, basis functions
integer :: e                     ! element index
integer :: i,j,k,l                   ! basis fn indices
integer :: al,be                 ! "alpha", "beta": local basis fn indices
integer :: iq                    ! quadrature point index
real(dp) :: jac                  ! Jacobian of transformation from parent
    ! coords xi in [-1,1] to coords x in [xa,xb]: x = (xb-xa)/2 xi + (xb+xa)/2
real(dp) :: phihq(size(xiq),size(xin)), dphihq(size(xiq),size(xin))   ! parent
    ! basis fns and derivs at quadrature points:
    ! phihq(i,j) = value of jth function at ith quadrature point
real(dp), dimension(size(xiq), size(ib,2)) :: f, xq
real(dp), dimension(size(xiq), size(ib,2), size(xin)*(size(xin)+1)/2, &
    size(ib,2), 0:max_k) :: Ykq
real(dp), dimension(size(xin)*(size(xin)+1)/2, size(ib,2), &
    size(xin)*(size(xin)+1)/2, size(ib,2), 0:max_k) :: mat2
integer :: e2
integer :: n, m, k_, nn

! initializations
Ne = size(ib, 2)
Nb = maxval(ib)
p = size(xin)-1
if (size(xin) /= size(ib,1)) &
    call stop_error("Error: inconsistent parent node and connectivity dimensions.")
! tabulate parent basis and derivatives at quadrature points
do al = 1, p+1
   do iq = 1, size(xiq)
      phihq(iq, al) = phih(xin, al, xiq(iq))
      dphihq(iq, al) = dphih(xin, al, xiq(iq))
   end do
end do
call get_quad_pts(xe, xiq, xq)

print *, "SIZE Ykq : ", str(size(Ykq ) * 8 * 1e-6_dp, 2), " MB"
print *, "SIZE mat2: ", str(size(mat2) * 8 * 1e-6_dp, 2), " MB"
print *, "Calculating Hartree screening functions..."
do e = 1, Ne
    print *, "    element =", e
    do j = 1, p+1
        do i = 1, j
            f = 0
            f(:, e) = phihq(:, i) * phihq(:, j)
            do k = 0, max_k
                Ykq(:, :, sym_idx(i, j), e, k) = hartree_y(k, p, xe, xiq, &
                    wtq, f) / xq
            end do
        end do
    end do
end do
print *, "Calculating 2 particle matrix elements..."
do e = 1, Ne
    print *, "    element =", e
    jac=(xe(e+1)-xe(e))/2
    do j = 1, p+1
        do i = 1, j
            do e2 = 1, Ne
                do n = 1, p+1
                    do m = 1, n
                        ! Note: we should use the fact that the elements
                        ! sym_idx(m, n), e2, sym_idx(i, j), e
                        ! sym_idx(i, j), e, sym_idx(m, n), e2
                        ! are equal (all other symmetries are already
                        ! exploited)
                        do k = 0, max_k
                            mat2(sym_idx(m, n), e2, sym_idx(i, j), e, k) = &
                                sum(jac * wtq * phihq(:, i) * phihq(:, j) *&
                                    Ykq(:, e, sym_idx(m, n), e2, k))
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
! accumulate Am matrix
! compute lower triangle
int2 = 0
print *, "Calculating non-local term..."
do e = 1, Ne
    print *, "    element =", e
    do al = 1, p+1
        i = ib(al, e)
        if (i == 0) cycle
        do nn = 1, p+1
            j = ib(nn, e)
            if (j == 0) cycle
            if (i < j) cycle
            do e2 = 1, Ne
                do be = 1, p+1
                    k = ib(be, e2)
                    if (k == 0) cycle
                    do m = 1, p+1
                        l = ib(m, e2)
                        if (l == 0) cycle
                        if (k < l) cycle
                        if ((i-1)*i/2 + j < (k-1)*k/2 + l) cycle
                        do k_ = 0, max_k
                            int2(ijkl2intindex(i, j, k, l), k_) &
                                = int2(ijkl2intindex(i, j, k, l), k_) &
                                + mat2(sym_idx(m, be), e2, &
                                    sym_idx(nn, al), e, k_)
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
print *, "  Done."
end subroutine

end module
