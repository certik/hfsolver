module schroed_util
use types, only: dp
use linalg, only: eigh, eigvals
use utils, only: stop_error
implicit none
private
public lho, radial

contains

subroutine lho(Nb, xq, wq, B, Bp)
integer, intent(in) :: Nb
real(dp), intent(in) :: xq(:), wq(:),  B(:,:), Bp(:,:)
real(dp), allocatable :: Am(:,:), Bm(:,:), c(:,:), lam(:), hq(:)
real(dp) :: En
integer :: i, j
allocate(hq(size(xq)))
allocate(Am(Nb,Nb), Bm(Nb,Nb), c(Nb,Nb), lam(Nb))
print *, "Assembly"
! Construct matrices A and B
do i = 1, Nb
    do j = 1, Nb
        ! A
        ! Both of these work:
        hq = Bp(:,i)*Bp(:,j)/2 + B(:,i)*B(:,j)*(xq**2)/2
        Am(i,j) = sum(wq*hq)

        ! B
        hq = B(:,i)*B(:,j)
        Bm(i,j) = sum(wq*hq)
    end do
end do

print *, "Checking symmetry"
do j = 1, Nb
    do i = 1, j-1
        if (max(abs(Am(i,j)), abs(Am(j,i))) > tiny(1._dp)) then
            if (abs(Am(i,j)-Am(j,i)) / max(abs(Am(i,j)), abs(Am(j,i))) &
                    > 1e-8_dp) then
                print *, i, j, Am(i,j)-Am(j,i), Am(i,j), Am(j,i)
                call stop_error("Am not symmetric")
            end if
        end if
        if (abs(Bm(i,j)-Bm(j,i)) > 1e-12_dp) call stop_error("Bm not symmetric")
   end do
end do


print *, "Eigensolver"

lam = eigvals(Am)
print "('cond A: ', es10.2)", maxval(abs(lam))/minval(abs(lam))
lam = eigvals(Bm)
print "('cond B: ', es10.2)", maxval(abs(lam))/minval(abs(lam))

! Solve an eigenproblem
call eigh(Am, Bm, lam, c)

print *, "n, energy, exact energy, error"
do i = 1, min(Nb, 20)
    En = 0.5_dp + (i-1)
    print "(i4, f30.8, f18.8, es12.2)", i, lam(i), En, abs(lam(i)-En)
end do
end subroutine

subroutine radial(Nb, xq, wq, B, Bp, Z, l)
integer, intent(in) :: Nb, Z, l
real(dp), intent(in) :: xq(:), wq(:),  B(:,:), Bp(:,:)
real(dp), allocatable :: Am(:,:), Bm(:,:), c(:,:), lam(:), hq(:)
real(dp) :: En
integer :: i, j
allocate(hq(size(xq)))
allocate(Am(Nb,Nb), Bm(Nb,Nb), c(Nb,Nb), lam(Nb))
print *, "Assembly"
! Construct matrices A and B
do i = 1, Nb
    do j = 1, Nb
        ! A
        ! Both of these work:
        hq = Bp(:,i)*Bp(:,j)/2 + B(:,i)*B(:,j)*(-Z/xq+l*(l+1)/(2*xq**2))
        Am(i,j) = sum(wq*hq)

        ! B
        hq = B(:,i)*B(:,j)
        Bm(i,j) = sum(wq*hq)
    end do
end do

print *, "Checking symmetry"
do j = 1, Nb
    do i = 1, j-1
        if (max(abs(Am(i,j)), abs(Am(j,i))) > tiny(1._dp)) then
            if (abs(Am(i,j)-Am(j,i)) / max(abs(Am(i,j)), abs(Am(j,i))) &
                    > 1e-8_dp) then
                print *, i, j, Am(i,j)-Am(j,i), Am(i,j), Am(j,i)
                call stop_error("Am not symmetric")
            end if
        end if
        if (abs(Bm(i,j)-Bm(j,i)) > 1e-12_dp) call stop_error("Bm not symmetric")
   end do
end do


print *, "Eigensolver"
! Solve an eigenproblem
call eigh(Am, Bm, lam, c)

print *, "n, energy, exact energy, error"
do i = 1, min(Nb, 20)
    En = -Z**2/(2._dp*(i+l)**2)
    print "(i4, f30.8, f18.8, es12.2)", i, lam(i), En, abs(lam(i)-En)
end do
end subroutine

end module
