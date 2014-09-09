program test_poisson_fft
use types, only: dp
use constants, only: pi, i_
use ofdft_fft, only: reciprocal_space_vectors, real2fourier, integralG, &
        radial_potential_fourier
use utils, only: assert, linspace
use ewald_sums, only: ewald_box
implicit none
! Exact value for a point charge
real(dp), parameter :: E_ewald_exact = -6.3839193288315013_dp
! Converged value for a Gaussian charge with alpha=5 (should be close, but not
! exactly equal to E_ewald_exact)
real(dp), parameter :: E_gauss5_conv = -6.2425476594199409_dp

! Two charges:
real(dp), parameter :: E_ewald_exact2  = -9.1591267925365365_dp


call test1()
!call test2()
call test3()
call test4()
call test5()
call test6()
call test7()

contains

subroutine test1()
real(dp), allocatable :: G(:, :, :, :), G2(:, :, :), ne(:, :, :)
complex(dp), allocatable :: neG(:, :, :), VeeG(:, :, :)
real(dp) :: L, Eee, alpha, Z, X(3, 1)
integer :: Ng

L = 2
Ng = 32

allocate(ne(Ng, Ng, Ng), neG(Ng, Ng, Ng), VeeG(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng))
call reciprocal_space_vectors(L, G, G2)

alpha = 5
X(:, 1) = L/2
Z = 3
call gaussian_charges([Z], X, L, alpha, ne)
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
real(dp) :: L, Eee, X(3, 2), alpha, Z
integer :: Ng

L = 2
Ng = 32

allocate(ne(Ng, Ng, Ng), neG(Ng, Ng, Ng), VeeG(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng))
call reciprocal_space_vectors(L, G, G2)

alpha = 5
X(:, 1) = [L/4, L/4, L/4]
X(:, 2) = 3*[L/4, L/4, L/4]
Z = 3
call gaussian_charges([Z, -Z], X, L, alpha, ne)
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

subroutine gaussian_charges(q, r, L, alpha, ne)
real(dp), intent(in) :: q(:), r(:, :), L, alpha
real(dp), intent(out) :: ne(:, :, :)
integer :: i, j, k, a, b, c, n
real(dp) :: x(3), r2
ne = 0
do i = 1, size(ne, 1)
do j = 1, size(ne, 2)
do k = 1, size(ne, 3)
    x = ([i, j, k]-1) * L / shape(ne) ! x is in [0, L)^3
    do n = 1, size(r, 2)
        do a = -1, 1
        do b = -1, 1
        do c = -1, 1
            r2 = sum((x-r(:, n)+[a, b, c]*L)**2)
            ne(i, j, k) = ne(i, j, k) + &
                q(n)*alpha**3/pi**(3._dp/2)*exp(-alpha**2*r2)
        end do
        end do
        end do
    end do
end do
end do
end do
end subroutine

subroutine test6()
real(dp), allocatable :: G(:, :, :, :), G2(:, :, :), ne(:, :, :)
complex(dp), allocatable :: neG(:, :, :), VeeG(:, :, :)
integer, parameter :: natom = 8
real(dp) :: L, Eee, X(3, natom), xred(3, natom), alpha, q(natom), E_madelung
integer :: Ng, i

L = 2
Ng = 32

allocate(ne(Ng, Ng, Ng), neG(Ng, Ng, Ng), VeeG(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng))
call reciprocal_space_vectors(L, G, G2)

alpha = 5
! Cl^-
xred(:, 1) = [0._dp, 0._dp, 0._dp]
xred(:, 2) = [1._dp/2, 1._dp/2, 0._dp]
xred(:, 3) = [1._dp/2, 0._dp, 1._dp/2]
xred(:, 4) = [0._dp, 1._dp/2, 1._dp/2]
! Na^+
xred(:, 5) = [1._dp/2, 1._dp/2, 1._dp/2]
xred(:, 6) = [1._dp/2, 0._dp, 0._dp]
xred(:, 7) = [0._dp, 1._dp/2, 0._dp]
xred(:, 8) = [0._dp, 0._dp, 1._dp/2]
q = [-1, -1, -1, -1, 1, 1, 1, 1]*1._dp
X = xred*L
call gaussian_charges(q, X, L, alpha, ne)
call real2fourier(ne, neG)
VeeG = 4*pi*neG / G2
VeeG(1, 1, 1) = 0

Eee = integralG(L, VeeG*conjg(neG)) / 2
do i = 1, natom
    Eee = Eee - q(i)**2 * alpha / sqrt(2*pi) ! subtract self-energy
end do

alpha = 1.74756459463318219064_dp ! Madelung constant for NaCl
E_madelung = -2*alpha/L

Eee = Eee / 4

print *, "Ng = ", Ng
print *, "Eee       =", Eee
print *, "Eee_exact =", E_madelung
print *, "error     =", abs(Eee - E_madelung)

call assert(abs(Eee - E_madelung) < 5e-6_dp)
end subroutine

subroutine test7()
real(dp), allocatable :: G(:, :, :, :), G2(:, :, :), ne(:, :, :), fac(:, :, :),&
    R(:)
real(dp), allocatable :: ne0(:, :, :), Vee0G(:, :, :)
complex(dp), allocatable :: neG(:, :, :), VeeG(:, :, :), ne0G(:, :, :)
integer, parameter :: natom = 8
real(dp) :: L, Eee, X(3, natom), xred(3, natom), alpha, q(natom), E_madelung
real(dp) :: stress(6), forces(3, natom), forces_ewald(3, natom), rel
integer :: Ng, i

L = 2
Ng = 64

allocate(ne(Ng, Ng, Ng), neG(Ng, Ng, Ng), VeeG(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng), fac(Ng, Ng, Ng))
allocate(ne0(Ng, Ng, Ng), ne0G(Ng, Ng, Ng), Vee0G(Ng, Ng, Ng))
call reciprocal_space_vectors(L, G, G2)

alpha = 14
! Cl^-
xred(:, 1) = [0._dp, 0._dp, 0._dp]
xred(:, 2) = [1._dp/2, 1._dp/2, 0._dp]
xred(:, 3) = [1._dp/2, 0._dp, 1._dp/2]
xred(:, 4) = [0._dp, 1._dp/2, 1._dp/2]
! Na^+
xred(:, 5) = [1._dp/2, 1._dp/2, 1._dp/2]
xred(:, 6) = [1._dp/2, 0._dp, 0._dp]
xred(:, 7) = [0._dp, 1._dp/2, 0._dp]
xred(:, 8) = [0._dp, 0._dp, 1._dp/2]

! We shift the y-position here, this makes the dipole moment non-zero
xred(2, 2) = 1._dp/4

q = [-1, -1, -1, -1, 1, 1, 1, 1]*1._dp
X = xred*L
call gaussian_charges(q, X, L, alpha, ne)
call real2fourier(ne, neG)
VeeG = 4*pi*neG / G2
VeeG(1, 1, 1) = 0

Eee = integralG(L, VeeG*conjg(neG)) / 2
do i = 1, natom
    Eee = Eee - q(i)**2 * alpha / sqrt(2*pi) ! subtract self-energy
end do

alpha = 1.8285774522233_dp ! Madelung constant
E_madelung = -2*alpha/L

Eee = Eee / 4

print *, "Ng = ", Ng
print *, "Eee       =", Eee
print *, "Eee_exact =", E_madelung
print *, "error     =", abs(Eee - E_madelung)

call assert(abs(Eee - E_madelung) < 2e-11_dp)

call ewald_box(L, X, q, Eee, forces_ewald, stress)
call assert(abs(Eee / 4 - E_madelung) < 1e-13)
print *, "Ewald Forces:"
do i = 1, natom
    print *, i, forces_ewald(:, i)
end do
print *, "FFT Forces:"
allocate(R(10000))
R = linspace(1._dp/10000, 0.1_dp, 10000)
call radial_potential_fourier(R, erf(alpha*R)/R, L, 1._dp, Vee0G)
do i = 1, natom
    ! Important: in the shift theorem, we need to use plus sign in exp(+iG*X),
    ! because our FFT uses minus sign in the forward real -> fourier transform.
    fac = q(i)*L**3*Vee0G*aimag(neG*exp(i_ * &
        (G(:,:,:,1)*X(1,i) + G(:,:,:,2)*X(2,i) + G(:,:,:,3)*X(3,i))))
    forces(1, i) = sum(G(:,:,:,1)*fac)
    forces(2, i) = sum(G(:,:,:,2)*fac)
    forces(3, i) = sum(G(:,:,:,3)*fac)
end do
do i = 1, natom
    print *, i, forces(:, i)
end do
print *, "Errors:"
do i = 1, natom
    rel = sqrt(sum((forces(:, i)-forces_ewald(:, i))**2) / &
            sum(forces_ewald(:, i)**2))
    print *, i, rel
    if (i == 2 .or. i == 6) then
        call assert(rel < 5e-6_dp)
    else
        call assert(rel < 1e-12_dp)
    end if
end do
end subroutine

end program
