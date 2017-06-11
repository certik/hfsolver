module radial_util
use types, only: dp
use linalg, only: eigh, eigvals
use utils, only: stop_error
implicit none
private
public lho, radial, radial_dirac

contains

subroutine lho(xq, wq, B, Bp, lam, condA, condB)
real(dp), intent(in) :: xq(:), wq(:),  B(:,:), Bp(:,:)
real(dp), allocatable, intent(out) :: lam(:)
real(dp), intent(out) :: condA, condB
real(dp), allocatable :: Am(:,:), Bm(:,:), c(:,:), hq(:)
real(dp) :: En
integer :: i, j, Nb
Nb = size(B, 2)
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
condA = maxval(abs(lam))/minval(abs(lam))
print "('cond A: ', es10.2)", condA
lam = eigvals(Bm)
condB = maxval(abs(lam))/minval(abs(lam))
print "('cond B: ', es10.2)", condB

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


subroutine radial_dirac(Nb, xq, wq, B, Bp, Z, kappa, c, solvePQ)
integer, intent(in) :: Nb, Z, kappa
real(dp), intent(in) :: xq(:), wq(:),  B(:,:), Bp(:,:), c
logical, intent(in) :: solvePQ ! .true. solve for P/Q, otherwise for g/f
real(dp), allocatable :: Am(:,:), Bm(:,:), sol(:,:), lam(:), hq(:)
real(dp), allocatable :: Bi(:), Bj(:), Bip(:), Bjp(:)
real(dp), allocatable :: Vq(:), Vqp(:)
real(dp) :: En
integer :: i, j, l, relat, Nq_total
Nq_total = size(xq)

allocate(Am(2*Nb,2*Nb), Bm(2*Nb,2*Nb), sol(2*Nb,2*Nb), lam(2*Nb))
allocate(hq(Nq_total))
allocate(Bi(Nq_total), Bj(Nq_total), Bip(Nq_total), Bjp(Nq_total), Vq(Nq_total))
allocate(Vqp(Nq_total))

print *, "Assembly"
! Construct matrices A and B
do i = 1, Nb
    do j = 1, Nb
        Bi = B(:,i)
        Bj = B(:,j)
        Bip = Bp(:,i)
        Bjp = Bp(:,j)
        Vq = -Z/xq
        Vqp = Z/xq**2
        ! A11
        if (solvePQ) then
            hq = c**2*Bip*Bjp+Bi*Bj*((Vq+c**2)**2+c**2*(kappa*(kappa+1)/xq**2))
        else
            hq = -c**2*(-Bip*xq**2*Bjp -kappa*(kappa+1)*Bi*Bj) &
                +c**4*Bi*xq**2*Bj+Bi*xq**2*Vq**2*Bj+2*c**2*Bi*xq**2*Vq*Bj
        end if
        Am(i,j) = sum(wq*hq)
        ! A12
        if (solvePQ) then
            hq = c*Vq*(+Bip*Bj-Bi*Bjp + 2*kappa/xq*Bi*Bj)
        else
            hq = -c*(2*Bi*xq**2*Vq*Bjp+2*(1-kappa)*Bi*xq*Vq*Bj+Bi*xq**2*Vqp*Bj)
        end if
        Am(i,j+Nb) = sum(wq*hq)
        ! A21
        if (solvePQ) then
            hq = c*Vq*(-Bip*Bj+Bi*Bjp + 2*kappa/xq*Bi*Bj)
        else
            hq = c*(2*Bi*xq**2*Vq*Bjp+2*(1+kappa)*Bi*xq*Vq*Bj+Bi*xq**2*Vqp*Bj)
        end if
        Am(i+Nb,j) = sum(wq*hq)
        ! A22
        if (solvePQ) then
            hq = c**2*Bip*Bjp+Bi*Bj*((Vq-c**2)**2+ &
                c**2*(-kappa*(-kappa+1)/xq**2))
        else
            hq = -c**2*(-Bip*xq**2*Bjp -(-kappa)*(-kappa+1)*Bi*Bj) &
                +c**4*Bi*xq**2*Bj+Bi*xq**2*Vq**2*Bj-2*c**2*Bi*xq**2*Vq*Bj
        end if
        Am(i+Nb,j+Nb) = sum(wq*hq)

        ! B11
        if (solvePQ) then
            hq = B(:,i)*B(:,j)
        else
            hq = B(:,i)*B(:,j)*xq**2
        end if
        Bm(i,j) = sum(wq*hq)
        ! B22
        Bm(i+Nb,j+Nb) = Bm(i,j)
    end do
end do

print *, "Checking symmetry"
do j = 1, 2*Nb
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
call eigh(Am, Bm, lam, sol)

print *, "n, energy, exact energy, error"
do i = 1, Nb
    if (kappa > 0) then
        l = kappa
        relat = 3
    else
        l = -kappa-1
        relat = 2
    end if
    En = E_nl(c, l+i, l, real(Z, dp), relat)
    lam(i) = sqrt(lam(i)) - c**2
    print "(i4, f30.8, f18.8, es12.2)", i, lam(i), En, abs(lam(i)-En)
end do

!do i = 1, Nb
!    solsP(:,i) = 0
!    do j = 1, Nb
!        solsP(:,i) = solsP(:,i) + sol(j,i)*B(:,j)
!    end do
!end do
!
!open(newunit=u, file="dirac_sol.txt", status="replace")
!write(u,*) xq
!do i = 1, Nb
!    write(u,*) solsP(:,i)
!end do
!close(u)

contains

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

end subroutine


end module
