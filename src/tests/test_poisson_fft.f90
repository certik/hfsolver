program test_poisson_fft
use types, only: dp
use constants, only: pi
use ofdft_fft, only: reciprocal_space_vectors, real2fourier, integralG
use utils, only: assert
implicit none

! All variables are in Hartree atomic units

real(dp), allocatable :: G(:, :, :, :), G2(:, :, :), ne(:, :, :)
complex(dp), allocatable :: neG(:, :, :), VeeG(:, :, :)
real(dp) :: L, Eee, x(3), alpha, r, charge_pos(3), Z, Eee_exact
integer :: Ng, i, j, k

L = 2
Ng = 32

allocate(ne(Ng, Ng, Ng), neG(Ng, Ng, Ng), VeeG(Ng, Ng, Ng))
allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng))
call reciprocal_space_vectors(L, G, G2)


alpha = 5
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
Eee_exact = 1.3010949954051885_dp

print *, "Ng = ", Ng
print *, "Eee       =", Eee
print *, "Eee_exact =", Eee_exact
print *, "error     =", abs(Eee - Eee_exact)

call assert(abs(Eee - Eee_exact) < 1e-12)

end program
