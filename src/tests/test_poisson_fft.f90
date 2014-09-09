program test_poisson_fft
use types, only: dp
use constants, only: pi
use ofdft_fft, only: reciprocal_space_vectors, real2fourier, integralG
use utils, only: assert
use ewald_sums, only: ewald_box
implicit none
! Exact value for a point charge
real(dp), parameter :: E_ewald_exact = -6.3839193288315013_dp
! Converged value for a Gaussian charge with alpha=5 (should be close, but not
! exactly equal to E_ewald_exact)
real(dp), parameter :: E_gauss5_conv = -6.2425476594175731_dp

! Two charges:
real(dp), parameter :: E_ewald_exact2  = -9.1591267925365365_dp


call test1()
!call test2()
call test3()
call test4()
call test5()

contains

subroutine test1()
real(dp), allocatable :: G(:, :, :, :), G2(:, :, :), ne(:, :, :)
complex(dp), allocatable :: neG(:, :, :), VeeG(:, :, :)
real(dp) :: L, Eee, x(3), alpha, r, charge_pos(3), Z
integer :: Ng, i, j, k

L = 2
Ng = 32

allocate(ne(Ng, Ng, Ng), neG(Ng, Ng, Ng), VeeG(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng))
call reciprocal_space_vectors(L, G, G2)

alpha = 5
Z = 3 ! Charge
charge_pos = L/2  ! The charge position is in the middle
do i = 1, size(ne, 1)
    do j = 1, size(ne, 2)
        do k = 1, size(ne, 3)
            x = ([i, j, k]-1) * L / shape(ne) ! x is in [0, L)^3
            r = sqrt(sum((x-charge_pos)**2))
            ne(i, j, k) = Z*alpha**3/pi**(3._dp/2)*exp(-alpha**2*r**2)
        end do
    end do
end do
call real2fourier(ne, neG)
VeeG = 4*pi*neG / G2
VeeG(1, 1, 1) = 0

Eee = integralG(L, VeeG*conjg(neG)) / 2
Eee = Eee - Z**2 * alpha / sqrt(2*pi) ! subtract self-energy

print *, "Ng = ", Ng
print *, "Eee       =", Eee
print *, "Eee_exact =", E_gauss5_conv
print *, "error     =", abs(Eee - E_gauss5_conv)

call assert(abs(Eee - E_gauss5_conv) < 2e-12_dp)
end subroutine


subroutine test2()
real(dp), allocatable :: G(:, :, :, :), G2(:, :, :), ne(:, :, :)
complex(dp), allocatable :: neG(:, :, :), VeeG(:, :, :)
real(dp) :: L, Eee, x(3), alpha, r, charge_pos(3), Z
integer :: Ng, i, j, k, ialpha

L = 2

print *, "Ng, alpha, Eee"
do Ng = 32, 256, 16
    allocate(ne(Ng, Ng, Ng), neG(Ng, Ng, Ng), VeeG(Ng, Ng, Ng))
    allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng))
    call reciprocal_space_vectors(L, G, G2)

    do ialpha = 1, 200, 10
        alpha = 5 + ialpha
        do i = 1, size(ne, 1)
            do j = 1, size(ne, 2)
                do k = 1, size(ne, 3)
                    x = ([i, j, k]-1) * L / shape(ne) ! x is in [0, L)^3
                    charge_pos = L/2  ! The charge position is in the middle
                    Z = 1 ! Charge
                    r = sqrt(sum((x-charge_pos)**2))
                    ne(i, j, k) = Z*alpha**3/pi**(3._dp/2)*exp(-alpha**2*r**2)
                end do
            end do
        end do
        call real2fourier(ne, neG)
        VeeG = 4*pi*neG / G2
        VeeG(1, 1, 1) = 0

        Eee = integralG(L, VeeG*conjg(neG)) / 2

        print *, Ng, alpha, Eee
    end do
    deallocate(ne, neG, VeeG, G, G2)
end do
end subroutine


subroutine test3()
integer :: N
real(dp) :: L, E_ewald, stress(6)
real(dp), allocatable :: X(:, :), fewald(:, :), q(:)
L = 2
N = 1
allocate(X(3, N), fewald(3, N), q(N))
X(:, 1) = [L/2, L/2, L/2]
q = 3
call ewald_box(L, X, q, E_ewald, fewald, stress)
print *, "Ewald       =", E_ewald
print *, "Ewald_exact =", E_ewald_exact
print *, "error       =", abs(E_ewald - E_ewald_exact)

call assert(abs(E_ewald - E_ewald_exact) < 1e-12_dp)
end subroutine

subroutine test4()
real(dp), allocatable :: G(:, :, :, :), G2(:, :, :), ne(:, :, :)
complex(dp), allocatable :: neG(:, :, :), VeeG(:, :, :)
real(dp) :: L, Eee, x(3), alpha, r, charge_pos(3), Z
integer :: Ng, i, j, k, a, b, c

L = 2
Ng = 32

allocate(ne(Ng, Ng, Ng), neG(Ng, Ng, Ng), VeeG(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng))
call reciprocal_space_vectors(L, G, G2)

alpha = 5
ne = 0
do i = 1, size(ne, 1)
    do j = 1, size(ne, 2)
        do k = 1, size(ne, 3)
            x = ([i, j, k]-1) * L / shape(ne) ! x is in [0, L)^3
            charge_pos = L/4
            Z = 3 ! Charge
            do a = -1, 1
            do b = -1, 1
            do c = -1, 1
                r = sqrt(sum((x-charge_pos+[a, b, c]*L)**2))
                ne(i, j, k) = ne(i, j, k) + &
                    Z*alpha**3/pi**(3._dp/2)*exp(-alpha**2*r**2)
            end do
            end do
            end do
            charge_pos = 3*L/4
            Z = -3 ! Charge
            do a = -1, 1
            do b = -1, 1
            do c = -1, 1
                r = sqrt(sum((x-charge_pos+[a, b, c]*L)**2))
                ne(i, j, k) = ne(i, j, k) + &
                    Z*alpha**3/pi**(3._dp/2)*exp(-alpha**2*r**2)
            end do
            end do
            end do
        end do
    end do
end do
call real2fourier(ne, neG)
VeeG = 4*pi*neG / G2
VeeG(1, 1, 1) = 0

Eee = integralG(L, VeeG*conjg(neG)) / 2
Eee = Eee - 2 * Z**2 * alpha / sqrt(2*pi) ! subtract self-energy (2 charges)

print *, "Ng = ", Ng
print *, "Eee       =", Eee
print *, "Eee_exact =", E_ewald_exact2
print *, "error     =", abs(Eee - E_ewald_exact2)

call assert(abs(Eee - E_ewald_exact2) < 1e-12_dp)
end subroutine

subroutine test5()
integer :: N
real(dp) :: L, E_ewald, stress(6)
real(dp), allocatable :: X(:, :), fewald(:, :), q(:)
L = 2
N = 2
allocate(X(3, N), fewald(3, N), q(N))
X(:, 1) = [L/4, L/4, L/4]
X(:, 2) = 3*[L/4, L/4, L/4]
q = [3, -3]
call ewald_box(L, X, q, E_ewald, fewald, stress)
print *, "Ewald       =", E_ewald
print *, "Ewald_exact =", E_ewald_exact2
print *, "error       =", abs(E_ewald - E_ewald_exact2)

call assert(abs(E_ewald - E_ewald_exact2) < 1e-12_dp)
end subroutine

end program
