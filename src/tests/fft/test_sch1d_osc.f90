program test_sch1d_osc

! Tests 1D harmonic oscillator ground state minimization (imaginary time
! propagation).

use types, only: dp
use constants, only: i_, pi
use ofdft, only: read_pseudo
use ofdft_fft, only: free_energy, radial_potential_fourier, &
    reciprocal_space_vectors, free_energy_min, real2fourier, integral, &
    fourier2real, real_space_vectors, vtk_save, integralG
use utils, only: loadtxt, stop_error, assert, linspace, strfmt
use splines, only: spline3pars, iixmin, poly3, spline3ders
use interp3d, only: trilinear
use md, only: positions_fcc, positions_bcc
use arpack, only: eig
implicit none
integer :: Ng
real(dp), allocatable :: G(:), G2(:)
real(dp), allocatable :: ne(:)
real(dp), allocatable :: Xn(:), Vn(:)
real(dp), allocatable :: d(:), v(:,:)
complex(dp), allocatable, dimension(:) :: psi, psiG
real(dp) :: L
integer :: i, j, nev, ncv
integer, parameter :: nelec = 1
real(dp) :: dt, psi_norm, E_tot
integer :: u, u2
real(dp) :: t
!real(dp), parameter :: Ng_list(*) = [2, 4, 8, 16, 64, 128, 256, 512]
real(dp), parameter :: Ng_list(*) = [2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, &
    20, 24, 25, 27, 30, 32, 36, 40, 45, 48, 50, 54, 60, 64, 72, 75, 80, 81, &
    90, 96, 100, 108, 120, 125, 128, 135, 144, 150, 160, 162, 180, 192, 200, &
    216, 225, 240, 243, 250, 256, 270, 288, 300, 320, 324, 360, 375, 384, 512, &
    1024]

L = 8


open(newunit=u2, file="pw.txt", status="replace")
do j = 1, size(Ng_list)
    Ng = Ng_list(j)
    allocate(ne(Ng))
    allocate(G(Ng), G2(Ng), psi(Ng))
    allocate(psiG(Ng))
    allocate(Xn(Ng), Vn(Ng))

    call real_space_vectors(L, Xn)
    call reciprocal_space_vectors(L, G, G2)

    open(newunit=u, file="sch1d_grid.txt", status="replace")
    write(u, *) Ng
    write(u, *) Xn
    Vn = gaussian_potential(Xn, 12._dp, L/2)
    write(u, *) Vn
    psi = gaussian_density(Xn, 12._dp, L/2)
    write(u, *) real(psi, dp)

    ! Solve Poisson
    call real2fourier(psi, psiG)
    psiG(1) = 0; psiG(2:) = 4*pi*psiG(2:) / G2(2:)
    call fourier2real(psiG, Vn)

    write(u, *) Vn
    close(u)

    nev = min(6, Ng-1)
    ncv = min(160, Ng)
    allocate(v(Ng,ncv), d(ncv))
    call eig(Ng, nev, ncv, "SA", av, d, v)
    print *, "n  eig  eig_integral"
    open(newunit=u, file="sch1d_psi.txt", status="replace")
    do i = 1, nev
        psi = v(:,i)
        ne = real(psi*conjg(psi), dp)
        psi_norm = integral(L, ne)
        psi = sqrt(nelec / psi_norm) * psi

        ne = real(psi*conjg(psi), dp)
        call real2fourier(psi, psiG)
        E_tot = 1._dp/2 * integralG(G2*abs(psiG)**2, L) + integral(L, Vn*ne)
        print *, i, d(i), E_tot
        write(u, *) real(psi, dp)
    end do
    close(u)

    write(u2,*) Ng, L, d(:nev)

    deallocate(ne, G, G2, psi, psiG, Xn, Vn, v, d)
end do
close(u)

contains

    pure function gaussian_density(x, alpha, x0) result(n)
    ! Returns a Gaussian charge.
    !
    ! This function returns the particle density `n` corresponding to the
    ! following radial potential `V`:
    !     V = -erf(alpha*r)/r
    ! by solving the radial Poisson equation (r*V(r))'' = -4*pi*r*n(r), e.g. in
    ! SymPy:
    !
    !     >>> var("alpha r", positive=True)
    !     >>> V = -erf(alpha*r)/r
    !     >>> n = (r*V).diff(r, 2) / (-4*pi*r)
    !     >>> n
    !     -alpha**3*exp(-alpha**2*r**2)/pi**(3/2)
    !
    ! Note that in 1D the potential `V` as calculated from the 1D Poisson
    ! equation V''(x) = -4*pi*n(x) is equal to:
    !
    !     >>> simplify(integrate(integrate(-4*pi*n, r), r))
    !     2*alpha**2*r*erf(alpha*r) + 2*alpha*exp(-alpha**2*r**2)/sqrt(pi)
    !
    ! The last potential is returned by gaussian_potential().
    real(dp), intent(in) :: x(:), alpha, x0
    real(dp) :: n(size(x))
    n = -alpha**3 * exp(-alpha**2*(x-x0)**2)/pi**(3._dp/2)
    end function

    pure function gaussian_potential(x, alpha, x0) result(V)
    ! Returns a Gaussian potential. See gaussian_density() for details.
    real(dp), intent(in) :: x(:), alpha, x0
    real(dp) :: V(size(x)), r(size(x))
    r = abs(x-x0)
    V = 2*alpha**2*r*erf(alpha*r) + 2*alpha*exp(-alpha**2*r**2)/sqrt(pi)
    end function

    subroutine av(x, y)
    ! Compute y = A*x
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    complex(dp), dimension(Ng) :: psi, psiG
    call real2fourier(x, psiG)
    call fourier2real(G2/2*psiG, psi)
    y = real(psi,dp) + Vn*x
    end

end program
