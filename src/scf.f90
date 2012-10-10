module scf
use types, only: dp
use linalg, only: eigh
use utils, only: stop_error, str
implicit none
private
public doscf, electronic_energy, get_nuclear_energy, ijkl2intindex, &
        ijkl2intindex2

contains

function formP(C, Noccupied) result(P)
! Forms the density matrix from occupied eigenvectors
real(dp), intent(in) :: C(:, :) ! C(i, j) is the i-th component
    ! of the j-th eigenvector
integer, intent(in) :: Noccupied ! the number of occupied orbitals
real(dp) :: P(size(C, 1), size(C, 1))
integer :: i, j
do i = 1, size(C, 1)
    do j = 1, size(C, 1)
        P(i, j) = 2 * sum(C(i, :Noccupied) * C(j, :Noccupied))
    end do
end do
end function

function formG(P, TT) result(G)
! Assembles Coulomb and exchange term
real(dp), intent(in) :: P(:, :) ! Density matrix
real(dp), intent(in) :: TT(:, :, :, :) ! T(k, l, i, j) = (ij|kl) - 0.5*(il|kj)
real(dp) :: G(size(P, 1), size(P, 1)) ! G(i,j)=sum_kl P(k,l) * T(k, l, i, j)
integer :: i, j
do i = 1, size(P, 1)
    do j = 1, size(P, 1)
        G(i, j) = sum(P * TT(:, :, i, j))
    end do
end do
end function

function formG2(P, int2) result(G)
! Assembles Coulomb and exchange term
! The same as formG, but doesn't require the big matrix TT(:, :, :, :)
real(dp), intent(in) :: P(:, :) ! Density matrix
real(dp), intent(in) :: int2(:)
real(dp), dimension(size(P, 1), size(P, 1)) :: G, TT
integer :: i, j, k, l, n
n = size(P, 1)
do i = 1, n
    do j = 1, n
        do k = 1, n
            do l = 1, n
                TT(k, l) = int2(ijkl2intindex(i, j, k, l)) &
                      - 0.5_dp * int2(ijkl2intindex(i, l, k, j))
            end do
        end do
        G(i, j) = sum(P * TT)
    end do
end do
end function

real(dp) function electronic_energy(P, H, F) result(E0)
! Calculates electronic energy
real(dp), intent(in) :: P(:, :), H(:, :), F(:, :)
E0 = 0.5_dp * sum(P*(H + F))
end function

real(dp) function kinetic_energy(P, T) result(Ekin)
! Calculates kinetic energy
real(dp), intent(in) :: P(:, :), T(:, :)
Ekin = sum(P * T)
end function

subroutine doscf(H, int2, S, Noccupied, Nscf, tolE, tolP, alpha, C, P, lam, E0)
! Runs the SCF cycle
real(dp), intent(in) :: H(:, :), int2(:), S(:, :)
integer, intent(in) :: Noccupied, Nscf
real(dp), intent(in) :: tolE, tolP, alpha
real(dp), intent(out) :: C(:, :), P(:, :), lam(:), E0
real(dp), dimension(size(H, 1), size(H, 1)) :: G, F, Pold
real(dp), dimension(size(H, 1), size(H, 1), size(H, 1), size(H, 1)) :: TT
integer :: i, j, k, l, n
real(dp) :: dP_, Escf(Nscf), Emin, Emax
n = size(H, 1)
print *, "Precalculating J-0.5*K ..."
do i = 1, n
    do j = 1, n
        do k = 1, n
            do l = 1, n
!                print *, i, j, k, l, int2(ijkl2intindex(i, j, k, l)), &
!                    int2(ijkl2intindex(i, l, k, j))
                TT(k, l, i, j) = int2(ijkl2intindex(i, j, k, l)) &
                      - 0.5_dp * int2(ijkl2intindex(i, l, k, j))
!                if (abs(int2(ijkl2intindex(i, j, k, l)) - &
!                    int2(ijkl2intindex(i, l, k, j))) > 1e-10) then
!                    print *, int2(ijkl2intindex(1, 1, 2, 2))
!                    print *, int2(ijkl2intindex(1, 2, 2, 1))
!                    print *, int2(ijkl2intindex(1, 2, 1, 2))
!                    print *, int2(ijkl2intindex(2, 1, 1, 2))
!                    print *, int2(ijkl2intindex(2, 1, 2, 1))
!                    print *, int2(ijkl2intindex(2, 2, 1, 1))
!                    stop "OK"
!                end if
            end do
        end do
    end do
end do
print *, "Done."
P = 0
do i = 1, Nscf
    G = formG(P, TT)
    !G = formG2(P, int2)
    F = H + G
    E0 = electronic_energy(P, H, F)
    Escf(i) = E0
    print *, i, E0
    call eigh(F, S, lam, C)
    Pold = P
    P = formP(C, Noccupied)
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

real(dp) function get_nuclear_energy(atno, xyz) result(Enuc)
integer, intent(in) :: atno(:)
real(dp), intent(in) :: xyz(:, :)
integer :: i, j
real(dp) :: R(3)
Enuc = 0
do i = 1, size(atno)
    do j = 1, i-1
        R = xyz(:, i) - xyz(:, j)
        Enuc = Enuc + atno(i)*atno(j) / sqrt(dot_product(R, R))
    end do
end do
end function

integer recursive function ijkl2intindex_slow(i, j, k, l) result(r)
! Index into the int2(:) array using the chemistry notation (ij|kl)
! Here is how to iterate over all (ij|kl) combinations:
!    ijkl = 1
!    do i = 1, n
!        do j = 1, i
!            do k = 1, n
!                do l = 1, k
!                    if ((i-1)*i/2 + j < (k-1)*k/2 + l) cycle
!                    call assert(ijkl2intindex(i, j, k, l) == ijkl)
!                    ijkl = ijkl + 1
!                end do
!            end do
!        end do
!    end do
integer, intent(in) :: i, j, k, l
integer :: ij, kl
if (i < j) then
    r = ijkl2intindex(j, i, k, l)
elseif (k < l) then
    r = ijkl2intindex(i, j, l, k)
else
    ij = (i-1)*i/2 + j
    kl = (k-1)*k/2 + l
    if (ij < kl) then
        r = (kl-1)*kl/2 + ij
    else
        r = (ij-1)*ij/2 + kl
    end if
end if
end function

integer pure function ijkl2intindex(i_, j_, k_, l_) result(r)
! Index into the int2(:) array using the chemistry notation (ij|kl)
! Here is how to iterate over all (ij|kl) combinations:
!    ijkl = 1
!    do i = 1, n
!        do j = 1, i
!            do k = 1, n
!                do l = 1, k
!                    if ((i-1)*i/2 + j < (k-1)*k/2 + l) cycle
!                    call assert(ijkl2intindex(i, j, k, l) == ijkl)
!                    ijkl = ijkl + 1
!                end do
!            end do
!        end do
!    end do
integer, intent(in) :: i_, j_, k_, l_
integer :: i, j, k, l
integer :: ij, kl
if (i_ < j_) then
    i = j_
    j = i_
else
    i = i_
    j = j_
end if
if (k_ < l_) then
    k = l_
    l = k_
else
    k = k_
    l = l_
end if
ij = (i-1)*i/2 + j
kl = (k-1)*k/2 + l
if (ij < kl) then
    r = (kl-1)*kl/2 + ij
else
    r = (ij-1)*ij/2 + kl
end if
end function

integer recursive function ijkl2intindex2(i, j, k, l) result(r)
! Index into the int2(:) array using the chemistry notation (ij|kl)
! Only the 4 general symmetries are assumed.
! Here is how to iterate over all (ij|kl) combinations:
!    ijkl = 1
!    do i = 1, n
!        do j = 1, i
!            do k = 1, i
!                do l = 1, i
!                    if (i == j .and. k < l) cycle
!                    if (i == k .and. j < l) cycle
!                    if (i == l .and. j < k) cycle
!                    call assert(ijkl2intindex2(i, j, k, l) == ijkl)
!                    ijkl = ijkl + 1
!                end do
!            end do
!        end do
!    end do
integer, intent(in) :: i, j, k, l
if (i < j .or. (i == j .and. k < l)) then
    r = ijkl2intindex2(j, i, l, k)
else if (i < k .or. (i == k .and. j < l)) then
    r = ijkl2intindex2(k, l, i, j)
else if (i < l .or. (i == l .and. j < k)) then
    r = ijkl2intindex2(l, k, j, i)
else
    r = getni(i) + getnj(i, j) + getnk(i, j, k) + getnl(i, j, k, l)
end if
end function

integer pure function getnl(i, j, k, l) result(n)
integer, intent(in) :: i, j, k, l
n = l
if (i == j) n = min(k, n)
if (i == k) n = min(j, n)
if (i <= n .and. j < k) n = n - 1
end function

integer pure function getnk(i, j, k_) result(n)
integer, intent(in) :: i, j, k_
integer :: k
n = 0
do k = 1, k_ - 1
    n = n + getnl(i, j, k, i)
end do
end function

integer pure function getnj(i, j_) result(n)
integer, intent(in) :: i, j_
integer :: j
n = 0
do j = 1, j_ - 1
    n = n + getnk(i, j, i+1)
end do
end function

integer pure function getni(i) result(n)
integer, intent(in) :: i
n = (i-1)**2 * ((i-1)**2 + 3) / 4
end function

end module
