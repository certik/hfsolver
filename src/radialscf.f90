module radialscf
use types, only: dp
use linalg, only: eigh, inv
use utils, only: stop_error, str, assert
use special_functions, only: wigner3j, getgaunt, getgauntr
use scf, only: ijkl2intindex, ijkl2intindex2
implicit none
private
public doscf, electronic_energy, kinetic_energy, slater2int2, get_basis, &
    radialC2C, radiallam2lam, slater2int22, numocc


! Array of wigner3j coefficients:
! wigner3j_(l1, l2, l3) = wigner3j(l1, l2, l3, 0, 0, 0)
! Use init_wigner3j() to initialize this array
real(dp), allocatable :: wigner3j_(:, :, :)

contains

integer function numocc(occ) result(k)
! returns index of last nonzero value in occ array
real(dp), intent(in) :: occ(:)
integer :: i
do i = size(occ), 1, -1
    if (abs(occ(i)) > 0) then
        k = i
        return
    end if
end do
k = 0
end function

function formP(nbfl, C, focc) result(P)
! Forms the density matrix from occupied eigenvectors
integer, intent(in) :: nbfl(0:) ! The length of each diagonal l-block
real(dp), intent(in) :: C(:, :, 0:) ! C(i, j) is the i-th component
    ! of the j-th eigenvector
real(dp), intent(in) :: focc(:, 0:)
real(dp) :: P(size(C, 1), size(C, 1), 0:ubound(nbfl, 1))
integer :: n, l, i, j
P = 0
do l = 0, ubound(nbfl, 1)
    do n = 1, numocc(focc(:, l))
        call assert(abs(focc(n, l) - 2*(2*l+1)) < tiny(1._dp))
        do j = 1, nbfl(l)
            do i = 1, j
                P(i, j, l) = P(i, j, l) + 2 * C(i, n, l) * C(j, n, l)
            end do
        end do
    end do

    ! Fill in the lower triangle:
    do j = 1, nbfl(l)
        do i = j+1, nbfl(l)
            P(i, j, l) = P(j, i, l)
        end do
    end do
end do
end function

function formG2(nbfl, P, slater) result(G)
! Assembles Coulomb and exchange term
! The same as formG, but doesn't require the big matrix TT(:, :, :, :)
integer, intent(in) :: nbfl(0:) ! The length of each diagonal l-block
real(dp), intent(in) :: P(:, :, 0:) ! Density matrix
real(dp), intent(in) :: slater(:, 0:)
real(dp), dimension(size(P, 1), size(P, 1), 0:ubound(nbfl, 1)) :: G
real(dp), dimension(size(P, 1), size(P, 1)) :: TT
integer :: l, lp, k, mu, nu, alpha, beta
real(dp) :: wig
G = 0
do l = 0, ubound(nbfl, 1)
  do nu = 1, nbfl(l)
    do mu = 1, nu
      do lp = 0, ubound(nbfl, 1)
        do beta = 1, nbfl(lp)
          do alpha = 1, nbfl(lp)
            TT(alpha, beta) = slater(ijkl2intindex(q(mu, l), q(nu, l), &
                q(beta, lp), q(alpha, lp)), 0)
            do k = abs(l-lp), abs(l+lp), 2
                wig = wigner3j_(l, k, lp)
                if (abs(wig) < tiny(1._dp)) cycle
                TT(alpha, beta) = TT(alpha, beta) - 0.5_dp * wig**2 &
                    * slater(ijkl2intindex(q(mu, l), q(alpha, lp), &
                        q(beta, lp), q(nu, l)), k)
            end do
          end do
        end do
        ! Note: The matrices G(mu, nu) and P(alpha, beta) are symmetric,
        ! but the matrix TT(alpha, beta) is *not* symmetric.
        G(mu, nu, l) = G(mu, nu, l) + (2*lp+1) * &
            sum(P(:nbfl(lp), :nbfl(lp), lp) * TT(:nbfl(lp), :nbfl(lp)))
      end do
    end do
  end do

  ! Fill in the lower triangle:
  do nu = 1, nbfl(l)
      do mu = nu+1, nbfl(l)
          G(mu, nu, l) = G(nu, mu, l)
      end do
  end do
end do

contains

integer function q(mu, l)
integer, intent(in) :: mu, l
q = sum(nbfl(:l-1)) + mu
end function

end function

function formG(nbfl, P, TT) result(G)
! Assembles Coulomb and exchange term
integer, intent(in) :: nbfl(0:) ! The length of each diagonal l-block
real(dp), intent(in) :: P(:, :, 0:) ! Density matrix
real(dp), intent(in) :: TT(:, :, 0:, :, :, 0:)
real(dp), dimension(size(P, 1), size(P, 1), 0:ubound(nbfl, 1)) :: G
integer :: l, lp, mu, nu
G = 0
do l = 0, ubound(nbfl, 1)
  do nu = 1, nbfl(l)
    do mu = 1, nu
      do lp = 0, ubound(nbfl, 1)
        G(mu, nu, l) = G(mu, nu, l) + (2*lp+1) * &
            sum(P(:nbfl(lp), :nbfl(lp), lp) * TT(:nbfl(lp), :nbfl(lp), &
            lp, mu, nu, l))
      end do
    end do
  end do

  ! Fill in the lower triangle:
  do nu = 1, nbfl(l)
      do mu = nu+1, nbfl(l)
          G(mu, nu, l) = G(nu, mu, l)
      end do
  end do
end do
end function

subroutine calcTT_indep(nbfl, slater, TT)
! Calculates the J-K term, only the upper triangle in (mu, nu) is calculated
integer, intent(in) :: nbfl(0:) ! The length of each diagonal l-block
real(dp), intent(in) :: slater(:, 0:)
! TT(mu, nu, l, alpha, beta, lp)
real(dp), intent(out) :: TT(:, :, 0:, :, :, 0:)
integer :: l, lp, k, mu, nu, alpha, beta
real(dp) :: wig
TT = 0
do l = 0, ubound(nbfl, 1)
  do nu = 1, nbfl(l)
    ! Important: only the upper triangle in (mu, nu) is calculated
    do mu = 1, nu
      do lp = 0, ubound(nbfl, 1)
        do beta = 1, nbfl(lp)
          do alpha = 1, nbfl(lp)
            TT(alpha, beta, lp, mu, nu, l) = &
                slater(ijkl2intindex(mu, nu, beta, alpha), 0)
            do k = abs(l-lp), abs(l+lp), 2
                wig = wigner3j_(l, k, lp)
                if (abs(wig) < tiny(1._dp)) cycle
                TT(alpha, beta, lp, mu, nu, l) = TT(alpha, beta, lp, mu, nu, l)&
                    - 0.5_dp * wig**2 &
                    * slater(ijkl2intindex(mu, alpha, beta, nu), k)
            end do
          end do
        end do
      end do
    end do
  end do
end do
end subroutine

function formG2_l_indep(nbfl, P, slater) result(G)
! Assembles Coulomb and exchange term
! The same as formG, but doesn't require the big matrix TT(:, :, :, :)
integer, intent(in) :: nbfl(0:) ! The length of each diagonal l-block
real(dp), intent(in) :: P(:, :, 0:) ! Density matrix
real(dp), intent(in) :: slater(:, 0:)
real(dp), dimension(size(P, 1), size(P, 1), 0:ubound(nbfl, 1)) :: G
real(dp), dimension(size(P, 1), size(P, 1)) :: TT
integer :: l, lp, k, mu, nu, alpha, beta
real(dp) :: wig
G = 0
do l = 0, ubound(nbfl, 1)
  do nu = 1, nbfl(l)
    do mu = 1, nu
      do lp = 0, ubound(nbfl, 1)
        do beta = 1, nbfl(lp)
          do alpha = 1, nbfl(lp)
            TT(alpha, beta) = slater(ijkl2intindex(mu, nu, beta, alpha), 0)
            do k = abs(l-lp), abs(l+lp), 2
                wig = wigner3j_(l, k, lp)
                if (abs(wig) < tiny(1._dp)) cycle
                TT(alpha, beta) = TT(alpha, beta) - 0.5_dp * wig**2 &
                    * slater(ijkl2intindex(mu, alpha, beta, nu), k)
            end do
          end do
        end do
        ! Note: The matrices G(mu, nu) and P(alpha, beta) are symmetric,
        ! but the matrix TT(alpha, beta) is *not* symmetric.
        G(mu, nu, l) = G(mu, nu, l) + (2*lp+1) * &
            sum(P(:nbfl(lp), :nbfl(lp), lp) * TT(:nbfl(lp), :nbfl(lp)))
      end do
    end do
  end do

  ! Fill in the lower triangle:
  do nu = 1, nbfl(l)
      do mu = nu+1, nbfl(l)
          G(mu, nu, l) = G(nu, mu, l)
      end do
  end do
end do
end function

real(dp) function electronic_energy(nbfl, P, H, F) result(E0)
! Calculates electronic energy
integer, intent(in) :: nbfl(0:) ! The length of each diagonal l-block
real(dp), intent(in) :: P(:, :, 0:), H(:, :, 0:), F(:, :, 0:)
integer :: l
E0 = 0
do l = 0, ubound(nbfl, 1)
    E0 = E0 + 0.5_dp * (2*l+1) * sum(P(:nbfl(l), :nbfl(l), l) &
         * (H(:nbfl(l), :nbfl(l), l) + F(:nbfl(l), :nbfl(l), l)))
end do
end function


real(dp) function kinetic_energy(nbfl, P, T) result(Ekin)
! Calculates kinetic energy
integer, intent(in) :: nbfl(0:) ! The length of each diagonal l-block
real(dp), intent(in) :: P(:, :, 0:), T(:, :, 0:)
integer :: l
Ekin = 0
do l = 0, ubound(nbfl, 1)
    Ekin = Ekin + (2*l+1) * sum(P(:nbfl(l), :nbfl(l), l) &
         * T(:nbfl(l), :nbfl(l), l))
end do
end function

subroutine doscf(nbfl, H, slater, S, focc, Nscf, tolE, tolP, alpha, &
        C, P, lam, E0, slater_l_indep, precalcTT, show_cond)
! Runs the SCF cycle
integer, intent(in) :: nbfl(0:) ! The length of each diagonal l-block
real(dp), intent(in) :: H(:, :, 0:), slater(:, 0:), S(:, :, 0:), focc(:, 0:)
integer, intent(in) :: Nscf
real(dp), intent(in) :: tolE, tolP, alpha
logical, intent(in), optional :: slater_l_indep ! the slater integrals are
! l-independent. Default: .false.
logical, intent(in), optional :: precalcTT ! The TT matrix will be
! precalculated (this needs much more memory, so for larger basis it might not
! be possible to use). Default: .false.
logical, intent(in), optional :: show_cond ! If present, then a condition
! number of the overlap matrix S will be shown for each "l".
real(dp), intent(out) :: C(:, :, 0:), P(:, :, 0:), lam(:, 0:), E0
real(dp), dimension(size(H, 1), size(H, 1), 0:ubound(nbfl, 1)) :: G, F, Pold
real(dp), allocatable :: TT(:, :, :, :, :, :)
integer :: i, l, n
real(dp) :: dP_, Escf(Nscf), Emin, Emax
logical :: l_indep, useTT
if (present(show_cond)) then
    print *, "Condition number of the overlap matrix:"
    print *, "l  cond"
    do l = 0, ubound(nbfl, 1)
        print "(i1, es10.2)", l, norm2(S(:nbfl(l), :nbfl(l), l)) * &
            norm2(inv(S(:nbfl(l), :nbfl(l), l)))
    end do
end if
l_indep = .false.
if (present(slater_l_indep)) l_indep = slater_l_indep
useTT = .false.
if (present(precalcTT)) useTT = precalcTT
call init_wigner3j(ubound(nbfl, 1), wigner3j_)
n = size(H, 1)
if (l_indep .and. useTT) then
    allocate(TT(n, n, 0:ubound(nbfl, 1), n, n, 0:ubound(nbfl, 1)))
    print *, "Calculating TT"
    call calcTT_indep(nbfl, slater, TT)
    print *, "Done"
end if
P = 0
do i = 1, Nscf
    if (l_indep) then
        if (useTT) then
            G = formG(nbfl, P, TT)
        else
            G = formG2_l_indep(nbfl, P, slater)
        end if
    else
        G = formG2(nbfl, P, slater)
    end if
    F = H + G
    E0 = electronic_energy(nbfl, P, H, F)
    Escf(i) = E0
    print *, i, E0
    lam = 0
    C = 0
    do l = 0, ubound(nbfl, 1)
        call eigh(F(:nbfl(l), :nbfl(l), l), S(:nbfl(l), :nbfl(l), l), &
            lam(:nbfl(l), l), C(:nbfl(l), :nbfl(l), l))
    end do
    Pold = P
    P = formP(nbfl, C, focc)
    dP_ = sqrt(sum((P - Pold)**2)/4)
    if (i > 3) then
        Emax = maxval(Escf(i-3:i))
        Emin = minval(Escf(i-3:i))
        if (dP_ < tolP .and. Emax-Emin < tolE) then
            write(*,*)  "Self-consistent convergence criterion achieved."
            print *, "SCF Iteration:", i
            print *, "tolE:", tolE
            print *, "tolP:", tolP
            print *, "dP:", dP_
            print *, "Emax-Emin:", Emax-Emin
            return
        end if
    end if
    ! Use linear mixing
    P = (1-alpha)*Pold + alpha*P
end do
call stop_error("SCF did not converge in " // str(Nscf) // " iterations.")
end subroutine

subroutine init_wigner3j(maxl, wigner3j_)
integer, intent(in) :: maxl
real(dp), intent(out), allocatable :: wigner3j_(:, :, :)
integer :: l, k, lp
allocate(wigner3j_(0:maxl, 0:2*maxl, 0:maxl))
do lp = 0, maxl
    do k = 0, 2*maxl
        do l = 0, maxl
            wigner3j_(l, k, lp) = wigner3j(l, k, lp, 0, 0, 0)
        end do
    end do
end do
end subroutine

subroutine slater2int2(nbfl, nlist, llist, mlist, slater, int2)
! Uses 8 symmetries. This is achieved by assuming the assumption 1):
!
! 1) The basis functions are real spherical harmonics
!
! With the assumption 1), the basis is real and thus the two particle integrals
! have 8 symmetries (just like for Cartesian gaussians). We then need integrals
! of 3 real spherical harmonics (so called real Gaunt coefficients).
!
! Unfortunately, this routine does not fully work, because the self energy then
! contains some spurious eigenvalues. Because the slater2int22 works correctly,
! this must be either due to a bug in real Gaunt coefficients (the getgauntr()
! routine is not well tested) or due to the fact that the assumption 1) is
! invalid --- this could happen because the rest of the code assumes complex
! spherical harmonics and thus there might be some mathematical inconsistency
! that renders this routine invalid.
integer, intent(in) :: nbfl(0:), nlist(:), llist(:), mlist(:)
real(dp), intent(in) :: slater(:, 0:)
real(dp), intent(out) :: int2(:)
integer :: a, b, c, d, Lmax, k_, q_
integer :: na, nb, nc, nd, la, lb, lc, ld, ma, mb, mc, md
integer :: ijkl
real(dp) :: r, fac
real(dp), allocatable :: ck(:, :, :, :, :)
real(dp), allocatable :: gr(:, :, :, :, :, :)
! Precalculate the radial integrals
Lmax = maxval(llist)
print *, "Calculating Gaunt coefficients..."
call getgaunt(Lmax, ck)
print *, "Calculating real Gaunt coefficients..."
call getgauntr(Lmax, gr)
!call print_gaunt(Lmax, ck)
print *, "Calculating int2..."
int2 = 0
ijkl = 0
do a = 1, size(nlist)
    print *, 100._dp * a / size(nlist), "%"
    do c = 1, a
        do b = 1, size(nlist)
            do d = 1, b
                if ((a-1)*a/2 + c < (b-1)*b/2 + d) cycle
                na = nlist(a); la = llist(a); ma = mlist(a)
                nb = nlist(b); lb = llist(b); mb = mlist(b)
                nc = nlist(c); lc = llist(c); mc = mlist(c)
                nd = nlist(d); ld = llist(d); md = mlist(d)
                ijkl = ijkl + 1
                r = 0
                do k_ = max(abs(la-lc), abs(lb-ld)), min(la+lc, lb+ld)
                    fac = 0
                    do q_ = -k_, k_
                        fac = fac + gr(la, ma, k_, q_, lc, mc) &
                                  * gr(lb, mb, k_, q_, ld, md)
                    end do
                    if (abs(fac) < tiny(1._dp)) cycle
                    r = r + fac * slater(ijkl2intindex(q(na, la), q(nc, lc), &
                            q(nb, lb), q(nd, ld)), k_)
                end do
                ! This must hold:
                !call assert(ijkl2intindex(a, c, b, d) == ijkl)
                int2(ijkl) = r
            end do
        end do
    end do
end do

contains

integer function q(mu, l)
integer, intent(in) :: mu, l
q = sum(nbfl(:l-1)) + mu
end function

end subroutine


subroutine slater2int22(nbfl, nlist, llist, mlist, slater, int2)
! Just like slater2int2, but only uses 4 symmetries. This version is well
! tested and it produces correct results.
integer, intent(in) :: nbfl(0:), nlist(:), llist(:), mlist(:)
real(dp), intent(in) :: slater(:, 0:)
real(dp), intent(out) :: int2(:)
integer :: a, b, c, d, Lmax, k_
integer :: na, nb, nc, nd, la, lb, lc, ld, ma, mb, mc, md
integer :: ijkl
real(dp) :: r, fac
real(dp), allocatable :: ck(:, :, :, :, :)
! Precalculate the radial integrals
Lmax = maxval(llist)
print *, "Calculating Gaunt coefficients..."
call getgaunt(Lmax, ck)
!call print_gaunt(Lmax, ck)
print *, "Calculating int2..."
int2 = 0
ijkl = 0
do a = 1, size(nlist)
    print *, 100._dp * a / size(nlist), "%"
    do c = 1, a
        do b = 1, a
            do d = 1, a
                if (a == c .and. b < d) cycle
                if (a == b .and. c < d) cycle
                if (a == d .and. c < b) cycle
                na = nlist(a); la = llist(a); ma = mlist(a)
                nb = nlist(b); lb = llist(b); mb = mlist(b)
                nc = nlist(c); lc = llist(c); mc = mlist(c)
                nd = nlist(d); ld = llist(d); md = mlist(d)
                ijkl = ijkl + 1
                if (ma+mb /= mc+md) cycle
                r = 0
                do k_ = max(abs(la-lc), abs(lb-ld)), min(la+lc, lb+ld)
                    fac = ck(k_, la, ma, lc, mc) * ck(k_, ld, md, lb, mb)
                    if (abs(fac) < tiny(1._dp)) cycle
                    r = r + fac * slater(ijkl2intindex(q(na, la), q(nc, lc), &
                            q(nb, lb), q(nd, ld)), k_)
                end do
                ! This must hold:
                call assert(ijkl2intindex2(a, c, b, d) == ijkl)
                int2(ijkl) = r
            end do
        end do
    end do
end do

contains

integer function q(mu, l)
integer, intent(in) :: mu, l
q = sum(nbfl(:l-1)) + mu
end function

end subroutine


subroutine slater2int2_4symmetries(nbfl, nlist, llist, mlist, slater, int2)
! This version uses complex spherical harmonics as a basis, thus
! the int2 integrals only have 4 symmetries. One should be using
! the ijkl2intindex2() here.
integer, intent(in) :: nbfl(0:), nlist(:), llist(:), mlist(:)
real(dp), intent(in) :: slater(:, 0:)
real(dp), intent(out) :: int2(:)
integer :: a, b, c, d, Lmax, k_
integer :: na, nb, nc, nd, la, lb, lc, ld, ma, mb, mc, md
integer :: ijkl
real(dp) :: r, fac
real(dp), allocatable :: ck(:, :, :, :, :)
! Precalculate the radial integrals
Lmax = maxval(llist)
print *, "Calculating Gaunt coefficients..."
call getgaunt(Lmax, ck)
!call print_gaunt(Lmax, ck)
print *, "Calculating int2..."
int2 = 0
ijkl = 0
do a = 1, size(nlist)
    print *, 100._dp * a / size(nlist), "%"
    do c = 1, a
        do b = 1, size(nlist)
            do d = 1, b
                if ((a-1)*a/2 + c < (b-1)*b/2 + d) cycle
                na = nlist(a); la = llist(a); ma = mlist(a)
                nb = nlist(b); lb = llist(b); mb = mlist(b)
                nc = nlist(c); lc = llist(c); mc = mlist(c)
                nd = nlist(d); ld = llist(d); md = mlist(d)
                ijkl = ijkl + 1
                if (ma+mb /= mc+md) cycle
                r = 0
                do k_ = max(abs(la-lc), abs(lb-ld)), min(la+lc, lb+ld)
                    fac = ck(k_, la, ma, lc, mc) * ck(k_, ld, md, lb, mb)
                    if (abs(fac) < tiny(1._dp)) cycle
                    r = r + fac * slater(ijkl2intindex(q(na, la), q(nc, lc), &
                            q(nb, lb), q(nd, ld)), k_)
                end do
                ! This must hold:
                !call assert(ijkl2intindex(a, c, b, d) == ijkl)
                int2(ijkl) = r
            end do
        end do
    end do
end do

contains

integer function q(mu, l)
integer, intent(in) :: mu, l
q = sum(nbfl(:l-1)) + mu
end function

end subroutine

subroutine get_basis(nbfl, nlist, llist, mlist)
integer, intent(in) :: nbfl(0:)
integer, allocatable, intent(out) :: nlist(:), llist(:), mlist(:)
integer :: i, n, l, m, a
a = 0
do l = 0, ubound(nbfl, 1)
    a = a + nbfl(l) * (2*l+1)
end do
allocate(nlist(a), llist(a), mlist(a))
i = 1
do l = 0, ubound(nbfl, 1)
    do n = 1, nbfl(l)
        do m = -l, l
            nlist(i) = n
            llist(i) = l
            mlist(i) = m
            i = i + 1
        end do
    end do
end do
end subroutine

subroutine radialC2C(nlist, llist, mlist, radialC, C)
integer, intent(in) :: nlist(:), llist(:), mlist(:)
real(dp), intent(in) :: radialC(:, :, 0:)
real(dp), intent(out) :: C(:, :)
integer :: a, i, na, la, ma, ni, li, mi
C = 0
do i = 1, size(nlist)
    ni = nlist(i); li = llist(i); mi = mlist(i)
    do a = 1, size(nlist)
        na = nlist(a); la = llist(a); ma = mlist(a)
        if (li == la .and. mi == ma) C(a, i) = radialC(na, ni, la)
    end do
end do
end subroutine

subroutine radiallam2lam(nlist, llist, radiallam, lam)
integer, intent(in) :: nlist(:), llist(:)
real(dp), intent(in) :: radiallam(:, 0:)
real(dp), intent(out) :: lam(:)
integer :: a, na, la
do a = 1, size(nlist)
    na = nlist(a); la = llist(a);
    lam(a) = radiallam(na, la)
end do
end subroutine

end module
