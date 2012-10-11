module mbpt
use types, only: dp
use scf, only: ijkl2intindex, ijkl2intindex2
use utils, only: assert
implicit none
private
public mbpt2, transform_int2, transform_int22, mbpt3, mbpt4

contains

function transform_int2(int2, C) result(moint2)
! Convert integrals 'int2' over gaussians into integrals over orbitals 'moint2'
! For that, just the orbitals 'C' are needed.
! int2(ijkl2intindex(mu, nu, alpha, beta)) is the (mu nu | alpha beta) integral
! over Gaussians:
real(dp), intent(in) :: int2(:)
! C(:, i) are the coefficients of the i-th orbital
real(dp), intent(in) :: C(:, :)
! int2(ijkl2intindex(i, j, k, l)) is the (i j | k l) integral over orbitals:
real(dp) :: moint2(size(int2))
integer :: mu, nu, alpha, beta ! Gaussian indices
integer :: i, j, k, l          ! orbital indices
integer :: ij, kl
integer :: n
real(dp) :: temp(size(C, 1), size(C, 1), size(C, 1), size(C, 1))
real(dp) :: temp2(size(C, 1), size(C, 1), size(C, 1), size(C, 1))
real(dp) :: tempvec(size(C, 1))
n = size(C, 1)
print *, "Transforming Gaussian (mu nu | alpha beta) integrals into"
print *, "orbital (i j | k l) integrals..."
print *, "  (mu nu | alpha beta) -> (mu nu | k beta)"
do mu = 1, n
    do nu = 1, n
        do beta = 1, n
            do alpha = 1, n
                tempvec(alpha) = int2(ijkl2intindex(mu, nu, alpha, beta))
            end do
            do k = 1, n
                temp(mu, nu, k, beta) = dot_product(C(:, k), tempvec)
            end do
        end do
    end do
end do
print *, "  (mu nu |   k   beta) -> ( i nu | k beta)"
do nu = 1, n
    do beta = 1, n
        do k = 1, n
            do i = 1, n
                temp2(i, nu, k, beta) = dot_product(C(:, i), temp(:, nu, k, beta))
            end do
        end do
    end do
end do
print *, "  ( i nu |   k   beta) -> ( i nu | k  l  )"
do i = 1, n
    do nu = 1, n
        do k = 1, n
            do l = 1, k  ! we can use kl symmetry here
                temp(i, nu, k, l) = dot_product(C(:, l), temp2(i, nu, k, :))
            end do
        end do
    end do
end do
print *, "  ( i nu |   k    l  ) -> ( i  j | k  l  )"
do k = 1, n
    do l = 1, k
        kl = (k-1)*k/2 + l
        do j = 1, n
            do i = j, n
                ij = (i-1)*i/2 + j
                if (ij < kl) cycle
                moint2(ijkl2intindex(i, j, k, l)) = &
                    dot_product(C(:, j), temp(i, :, k, l))
            end do
        end do
    end do
end do
print *, "Done."
end function

function transform_int22(int2, ijkl2intindex, C) result(moint2)
! Convert integrals 'int2' over gaussians into integrals over orbitals 'moint2'
! For that, just the orbitals 'C' are needed.
! int2(ijkl2intindex2(mu, nu, alpha, beta)) is the (mu nu | alpha beta) integral
! over Gaussians. Only 4 symmetries are assumed.
real(dp), intent(in) :: int2(:)
integer, intent(in) :: ijkl2intindex(:, :, :, :)
! C(:, i) are the coefficients of the i-th orbital
real(dp), intent(in) :: C(:, :)
! int2(ijkl2intindex(i, j, k, l)) is the (i j | k l) integral over orbitals:
real(dp) :: moint2(size(int2))
integer :: mu, nu, alpha, beta ! Gaussian indices
integer :: i, j, k, l          ! orbital indices
integer :: n, ijkl
real(dp) :: temp(size(C, 1), size(C, 1), size(C, 1), size(C, 1))
real(dp) :: temp2(size(C, 1), size(C, 1), size(C, 1), size(C, 1))
real(dp) :: tempvec(size(C, 1))
n = size(C, 1)
print *, "Transforming Gaussian (mu nu | alpha beta) integrals into"
print *, "orbital (i j | k l) integrals..."
print *, "  (mu nu | alpha beta) -> (mu nu | k beta)"
do mu = 1, n
    do nu = 1, n
        do beta = 1, n
            do alpha = 1, n
                tempvec(alpha) = int2(ijkl2intindex(mu, nu, alpha, beta))
            end do
            do k = 1, n
                temp(mu, nu, k, beta) = dot_product(C(:, k), tempvec)
            end do
        end do
    end do
end do
print *, "  (mu nu |   k   beta) -> ( i nu | k beta)"
do nu = 1, n
    do beta = 1, n
        do k = 1, n
            do i = 1, n
                temp2(beta, i, nu, k) = dot_product(C(:, i), temp(:, nu, k, beta))
            end do
        end do
    end do
end do
print *, "  ( i nu |   k   beta) -> ( i nu | k  l  )"
do i = 1, n
    do nu = 1, n
        do k = 1, n
            tempvec = temp2(:, i, nu, k)
            do l = 1, n
                temp(nu, k, l, i) = dot_product(C(:, l), tempvec)
            end do
        end do
    end do
end do
print *, "  ( i nu |   k    l  ) -> ( i  j | k  l  )"
ijkl = 1
do i = 1, n
    do j = 1, i
        do k = 1, i
            do l = 1, i
                if (i == j .and. k < l) cycle
                if (i == k .and. j < l) cycle
                if (i == l .and. j < k) cycle
                ! This must hold:
                !call assert(ijkl2intindex(i, j, k, l) == ijkl)
                moint2(ijkl) = dot_product(C(:, j), temp(:, k, l, i))
                ijkl = ijkl + 1
            end do
        end do
    end do
end do
print *, "Done."
end function

real(dp) function mbpt2(moint2, lam, Noccupied) result(E2)
use mbpt2, only: term1
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
E2 = term1(moint2, lam, Noccupied)
end function

real(dp) function mbpt3(moint2, lam, Noccupied) result(E3)
use mbpt3, only: term1, term2, term3
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
E3 = &
    + term1(moint2, lam, Noccupied) &
    + term2(moint2, lam, Noccupied) &
    + term3(moint2, lam, Noccupied)
end function

real(dp) function mbpt4(moint2, lam, Noccupied) result(E3)
use mbpt4
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
E3 =      term1(moint2, lam, Noccupied)
E3 = E3 + term2(moint2, lam, Noccupied)
E3 = E3 + term3(moint2, lam, Noccupied)
E3 = E3 + term4(moint2, lam, Noccupied)
E3 = E3 + term5(moint2, lam, Noccupied)
E3 = E3 + term6(moint2, lam, Noccupied)
E3 = E3 + term7(moint2, lam, Noccupied)
E3 = E3 + term8(moint2, lam, Noccupied)
E3 = E3 + term9(moint2, lam, Noccupied)
E3 = E3 + term10(moint2, lam, Noccupied)
E3 = E3 + term11(moint2, lam, Noccupied)
E3 = E3 + term12(moint2, lam, Noccupied)
E3 = E3 + term13(moint2, lam, Noccupied)
E3 = E3 + term14(moint2, lam, Noccupied)
E3 = E3 + term15(moint2, lam, Noccupied)
E3 = E3 + term16(moint2, lam, Noccupied)
E3 = E3 + term17(moint2, lam, Noccupied)
E3 = E3 + term18(moint2, lam, Noccupied)
E3 = E3 + term19(moint2, lam, Noccupied)
E3 = E3 + term20(moint2, lam, Noccupied)
E3 = E3 + term21(moint2, lam, Noccupied)
E3 = E3 + term22(moint2, lam, Noccupied)
E3 = E3 + term23(moint2, lam, Noccupied)
E3 = E3 + term24(moint2, lam, Noccupied)
E3 = E3 + term25(moint2, lam, Noccupied)
E3 = E3 + term26(moint2, lam, Noccupied)
E3 = E3 + term27(moint2, lam, Noccupied)
E3 = E3 + term28(moint2, lam, Noccupied)
E3 = E3 + term29(moint2, lam, Noccupied)
E3 = E3 + term30(moint2, lam, Noccupied)
E3 = E3 + term31(moint2, lam, Noccupied)
E3 = E3 + term32(moint2, lam, Noccupied)
E3 = E3 + term33(moint2, lam, Noccupied)
E3 = E3 + term34(moint2, lam, Noccupied)
E3 = E3 + term35(moint2, lam, Noccupied)
E3 = E3 + term36(moint2, lam, Noccupied)
E3 = E3 + term37(moint2, lam, Noccupied)
E3 = E3 + term38(moint2, lam, Noccupied)
E3 = E3 + term39(moint2, lam, Noccupied)
end function

end module
