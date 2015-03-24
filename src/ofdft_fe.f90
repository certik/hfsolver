module ofdft_fe

use types, only: dp
use feutils, only: phih, dphih
use fe_mesh, only: cartesian_mesh_3d, define_connect_tensor_3d, &
    c2fullc_3d, fe2quad_3d, vtk_save, fe_eval_xyz, line_save, &
    fe2quad_3d_lobatto, quad2fe_3d
use poisson3d_assembly, only: assemble_3d, integral, func2quad, func_xyz, &
    assemble_3d_precalc, assemble_3d_csr, assemble_3d_coo_A, &
    assemble_3d_coo_rhs, local_overlap_matrix, assemble_3d_coo_rhs_spectral
use feutils, only: get_parent_nodes, get_parent_quad_pts_wts, quad_gauss, &
    quad_lobatto
!use linalg, only: solve
use ofdft, only: f
use isolve, only: solve_cg
use utils, only: assert, zeros, stop_error
use sparse, only: coo2csr_canonical
use constants, only: pi
use xc, only: xc_pz
use optimize, only: bracket, brent, parabola_vertex
use umfpack, only: factorize, solve, free_data, umfpack_numeric
use integration, only: integrate_trapz_1
implicit none
private
public free_energy_min, radial_density_fourier, fe_data, initialize_fe, &
    free_energy_min_low_level, WITH_UMFPACK, free_energy, free_energy2, &
    free_energy2_low_level

logical, parameter :: WITH_UMFPACK=.false.

type fe_data
    real(dp) :: Lx, Ly, Lz, jac_det
    integer :: p, Nb, Nq, Ne
    integer :: Nx, Ny, Nz ! Number of elements in each direction
    ! spetral==.true. if we using spectral elements (happens if and only if
    ! Nq=p+1 and quad_type==quad_lobatto)
    logical :: spectral
    ! nodes(:, i) are the (x,y) coordinates of the i-th mesh node
    real(dp), allocatable :: nodes(:, :)
    ! elems(:, i) are nodes of the i-th element
    integer, allocatable :: elems(:, :)
    real(dp), allocatable :: xin(:), xiq(:), wtq3(:, :, :), phihq(:, :)
    real(dp), allocatable :: phi_v(:, :, :, :, :, :), dphihq(:, :)
    integer, allocatable :: in(:, :, :, :), ib(:, :, :, :)
    ! The CSR matrix S (Sp, Sj, Sx) is the overlap matrix
    integer, allocatable :: Sp(:), Sj(:)
    real(dp), allocatable :: Sx(:)
    ! The CSR matrix A (Ap, Aj, Ax) is the dicrete Poisson system matrix
    integer, allocatable :: Ap(:), Aj(:)
    real(dp), allocatable :: Ax(:)
    ! If WITH_UMFPACK==.true., then 'matd' contains the factorized matrix A
    ! (Ap, Aj, Ax). Otherwise it is unused.
    type(umfpack_numeric) :: matd
end type

contains

subroutine free_energy_min(Nelec, L, Nex, Ney, Nez, p, T_au, fnen, fn_pos, &
        Nq, quad_type, energy_eps, &
        Eh, Een, Ts, Exc, Nb)
! Higher level subroutine that does FE precalculation inside it, so it is slow,
! but the interface is simpler.
real(dp), intent(in) :: Nelec, L, T_au
integer, intent(in) :: p, Nq, quad_type, Nex, Ney, Nez
procedure(func_xyz) :: fnen ! (negative) ionic particle density
procedure(func_xyz) :: fn_pos ! (positive) electronic particle density
real(dp), intent(in) :: energy_eps
real(dp), intent(out) :: Eh, Een, Ts, Exc
integer, intent(out) :: Nb

real(dp), allocatable, dimension(:, :, :, :) :: nenq_pos, nq_pos
type(fe_data) :: fed

call initialize_fe(L, Nex, Ney, Nez, p, Nq, quad_type, fed)

allocate(nenq_pos(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(nq_pos(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
nenq_pos = func2quad(fed%nodes, fed%elems, fed%xiq, fnen)
nq_pos = func2quad(fed%nodes, fed%elems, fed%xiq, fn_pos)

call free_energy_min_low_level(Nelec, T_au, nenq_pos, nq_pos, energy_eps, &
        fed, Eh, Een, Ts, Exc)

Nb = fed%Nb
end subroutine

subroutine initialize_fe(L, Nex, Ney, Nez, p, Nq, quad_type, fed)
! Initializes the FE datastructures (returned in 'fed').
! The 'fed' derived type can be assumed to be read-only after this function
! returns it.
real(dp), intent(in) :: L ! Length of the box, assumed to be [0, L]^3
integer, intent(in) :: Nex, Ney, Nez ! # of finite elements in each direction
integer, intent(in) :: p ! polynomial order
! Quadrature order (1..64) and quadrature type (quad_lobatto or quad_gauss)
integer, intent(in) :: Nq, quad_type
type(fe_data), intent(out) :: fed ! FE datastructures

integer :: Nn, ibc
integer :: i, j, k
integer :: Ncoo
integer, allocatable :: matAi_coo(:), matAj_coo(:)
real(dp), allocatable :: matAx_coo(:), wtq(:), Am_loc(:, :, :, :, :, :)
integer :: idx

ibc = 3 ! Periodic boundary condition

fed%Lx = L
fed%Ly = L
fed%Lz = L

call cartesian_mesh_3d(Nex, Ney, Nez, &
    [0, 0, 0]*1._dp, [fed%Lx, fed%Ly, fed%Lz], fed%nodes, fed%elems)
fed%Nx = Nex
fed%Ny = Ney
fed%Nz = Nez
Nn = size(fed%nodes, 2)
fed%Ne = size(fed%elems, 2)
fed%Nq = Nq
fed%p = p
fed%spectral = (Nq == p+1 .and. quad_type == quad_lobatto)

print *, "Number of nodes:", Nn
print *, "Number of elements:", fed%Ne
print *, "Nq =", Nq
print *, "p =", p
if (fed%spectral) then
    print *, "Spectral elements: ON"
else
    print *, "Spectral elements: OFF"
end if

allocate(fed%xin(p+1))
call get_parent_nodes(2, p, fed%xin)
allocate(fed%xiq(Nq), wtq(Nq), fed%wtq3(Nq, Nq, Nq))
call get_parent_quad_pts_wts(quad_type, Nq, fed%xiq, wtq)
forall(i=1:Nq, j=1:Nq, k=1:Nq) fed%wtq3(i, j, k) = wtq(i)*wtq(j)*wtq(k)
allocate(fed%phihq(size(fed%xiq), size(fed%xin)))
allocate(fed%dphihq(size(fed%xiq), size(fed%xin)))
! Tabulate parent basis at quadrature points
forall(i=1:size(fed%xiq), j=1:size(fed%xin)) &
        fed%phihq(i, j) =  phih(fed%xin, j, fed%xiq(i))
forall(i=1:size(fed%xiq), j=1:size(fed%xin)) &
        fed%dphihq(i, j) = dphih(fed%xin, j, fed%xiq(i))
allocate(Am_loc(Nq, Nq, Nq, Nq, Nq, Nq))
allocate(fed%phi_v(Nq, Nq, Nq, p+1, p+1, p+1))
print *, "Precalculating element matrix"
call assemble_3d_precalc(p, Nq, &
    fed%nodes(1, fed%elems(7, 1)) - fed%nodes(1, fed%elems(1, 1)), &
    fed%nodes(2, fed%elems(7, 1)) - fed%nodes(2, fed%elems(1, 1)), &
    fed%nodes(3, fed%elems(7, 1)) - fed%nodes(3, fed%elems(1, 1)), &
    fed%wtq3, fed%phihq, fed%dphihq, fed%jac_det, Am_loc, fed%phi_v)

print *, "local to global mapping"
call define_connect_tensor_3d(Nex, Ney, Nez, p, 1, fed%in)
call define_connect_tensor_3d(Nex, Ney, Nez, p, ibc, fed%ib)
fed%Nb = maxval(fed%ib)
Ncoo = fed%Ne*(p+1)**6
allocate(matAi_coo(Ncoo), matAj_coo(Ncoo), matAx_coo(Ncoo))
print *, "Assembling matrix A"
call assemble_3d_coo_A(fed%Ne, fed%p, fed%ib, Am_loc, matAi_coo, matAj_coo, matAx_coo, idx)
print *, "COO -> CSR"
call coo2csr_canonical(matAi_coo(:idx), matAj_coo(:idx), matAx_coo(:idx), &
    fed%Ap, fed%Aj, fed%Ax)
print *, "CSR Matrix:"
print *, "    dimension:", fed%Nb
print *, "    number of nonzeros:", size(fed%Ax)
print "('     density:',f11.8,'%')", size(fed%Ax)*100._dp / real(fed%Nb, dp)**2
if (WITH_UMFPACK) then
    print *, "umfpack factorize"
    call factorize(fed%Nb, fed%Ap, fed%Aj, fed%Ax, fed%matd)
    print *, "done"
end if

!print *, "Precalculating overlap element matrix"
!call local_overlap_matrix(fed%wtq3, fed%jac_det, fed%phi_v, Am_loc)
!print *, "Assembling matrix S"
!call assemble_3d_coo_A(fed%Ne, fed%p, fed%ib, Am_loc, matAi_coo, matAj_coo, matAx_coo, idx)
!print *, "COO -> CSR"
!call coo2csr_canonical(matAi_coo(:idx), matAj_coo(:idx), matAx_coo(:idx), &
!    fed%Sp, fed%Sj, fed%Sx)
!print *, "CSR Matrix:"
!print *, "    dimension:", fed%Nb
!print *, "    number of nonzeros:", size(fed%Sx)
!print "('     density:',f11.8,'%')", size(fed%Sx)*100._dp / real(fed%Nb, dp)**2
end subroutine

subroutine free_fe_umfpack(fed)
type(fe_data), intent(inout) :: fed
if (WITH_UMFPACK) then
    call free_data(fed%matd)
end if
end subroutine

subroutine free_energy_min_low_level(Nelec, T_au, nenq_pos, nq_pos, energy_eps,&
        ! Geometry, finite element specification and internal datastructures
        fed, &
        ! Output arguments
        Eh, Een, Ts, Exc)
! Lower level subroutine that doesn't do any FE precalculation, so it is fast.
real(dp), intent(in) :: Nelec
real(dp), intent(in) :: nenq_pos(:, :, :, :)
real(dp), intent(inout) :: nq_pos(:, :, :, :)
real(dp), intent(in) :: T_au
real(dp), intent(in) :: energy_eps
type(fe_data), intent(in) :: fed
real(dp), intent(out) :: Eh, Een, Ts, Exc

real(dp), allocatable :: free_energies(:)
real(dp), allocatable, dimension(:, :, :, :) :: Hpsi, &
    psi, psi_, psi_prev, ksi, ksi_prev, phi, phi_prime, eta, Venq, &
    nenq_neutral
integer :: iter, max_iter
real(dp) :: mu, last3, brent_eps, free_energy_, &
    gamma_d, gamma_n, theta, theta_a, theta_b, theta_c, fa, fb, fc
real(dp) :: f2
real(dp) :: psi_norm
real(dp) :: background
real(dp), allocatable :: rhs(:), sol(:), fullsol(:), fullsol2(:)

brent_eps = 1e-3_dp
max_iter = 200

allocate(nenq_neutral(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(Venq(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(Hpsi(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(psi(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(psi_(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(psi_prev(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(phi(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(phi_prime(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(ksi(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(ksi_prev(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(eta(fed%Nq, fed%Nq, fed%Nq, fed%Ne))


! Calculate Venq
allocate(rhs(fed%Nb), sol(fed%Nb), fullsol(maxval(fed%in)))
allocate(fullsol2(maxval(fed%in)))
background = integral(fed%nodes, fed%elems, fed%wtq3, nenq_pos) / &
    (fed%Lx*fed%Ly*fed%Lz)
print *, "Total (negative) ionic charge: ", background * (fed%Lx*fed%Ly*fed%Lz)
print *, "Subtracting constant background (Q/V): ", background
nenq_neutral = nenq_pos - background
print *, "Assembling RHS..."
if (fed%spectral) then
    call assemble_3d_coo_rhs_spectral(fed%Ne, fed%p, 4*pi*nenq_neutral, &
    fed%jac_det, fed%wtq3, fed%ib, rhs)
else
    call assemble_3d_coo_rhs(fed%Ne, fed%p, 4*pi*nenq_neutral, fed%jac_det, &
        fed%wtq3, fed%ib, fed%phi_v, rhs)
end if
print *, "sum(rhs):    ", sum(rhs)
print *, "integral rhs:", integral(fed%nodes, fed%elems, fed%wtq3, nenq_neutral)
print *, "Solving..."
if (WITH_UMFPACK) then
    call solve(fed%Ap, fed%Aj, fed%Ax, sol, rhs, fed%matd)
else
    sol = solve_cg(fed%Ap, fed%Aj, fed%Ax, rhs, zeros(size(rhs)), 1e-12_dp, 400)
end if
print *, "Converting..."
call c2fullc_3d(fed%in, fed%ib, sol, fullsol)
if (fed%spectral) then
    call fe2quad_3d_lobatto(fed%elems, fed%xiq, fed%in, fullsol, Venq)
else
    call fe2quad_3d(fed%elems, fed%xin, fed%xiq, fed%phihq, fed%in, fullsol, Venq)
end if
print *, "Done"


psi = sqrt(nq_pos)
psi_norm = integral(fed%nodes, fed%elems, fed%wtq3, psi**2)
print *, "Initial norm of psi:", psi_norm
psi = sqrt(Nelec / psi_norm) * psi
psi_norm = integral(fed%nodes, fed%elems, fed%wtq3, psi**2)
print *, "norm of psi:", psi_norm
! This returns H[n] = delta F / delta n, we save it to the Hpsi variable to
! save space:
call free_energy(fed%nodes, fed%elems, fed%in, fed%ib, fed%Nb, fed%Lx, fed%Ly, fed%Lz, fed%xin, fed%xiq, fed%wtq3, T_au, &
    Venq, psi**2, fed%phihq, fed%Ap, fed%Aj, fed%Ax, fed%matd, fed%spectral, &
    fed%phi_v, fed%jac_det, &
    Eh, Een, Ts, Exc, free_energy_, Hpsi=Hpsi)
! Hpsi = H[psi] = delta F / delta psi = 2*H[n]*psi, due to d/dpsi = 2 psi d/dn
Hpsi = Hpsi * 2*psi
mu = 1._dp / Nelec * integral(fed%nodes, fed%elems, fed%wtq3, 0.5_dp * psi * Hpsi)
ksi = 2*mu*psi - Hpsi
phi = ksi
phi_prime = phi - 1._dp / Nelec *  integral(fed%nodes, fed%elems, fed%wtq3, phi * psi) * psi
eta = sqrt(Nelec / integral(fed%nodes, fed%elems, fed%wtq3, phi_prime**2)) * phi_prime
theta = pi/2
!print *, "Summary of energies [a.u.]:"
!print "('    Ts   = ', f14.8)", Ts
!print "('    Een  = ', f14.8)", Een
!print "('    Eee  = ', f14.8)", Eh
!print "('    Exc  = ', f14.8)", Exc
!print *, "   ---------------------"
!print "('    Etot = ', f14.8, ' a.u.')", free_energy_
allocate(free_energies(max_iter))
gamma_n = 0
do iter = 1, max_iter
    theta_a = 0
    theta_b = mod(theta, 2*pi)
    call bracket(func, theta_a, theta_b, theta_c, fa, fb, fc, 100._dp, 20, verbose=.false.)
    if (iter < 2) then
        call brent(func, theta_a, theta_b, theta_c, brent_eps, 50, theta, &
            free_energy_, verbose=.false.)
    else
        call parabola_vertex(theta_a, fa, theta_b, fb, theta_c, fc, theta, f2)
    end if
    ! TODO: We probably don't need to recalculate free_energy_ here:
    psi_prev = psi
    psi = cos(theta) * psi + sin(theta) * eta
    call free_energy(fed%nodes, fed%elems, fed%in, fed%ib, fed%Nb, fed%Lx, &
        fed%Ly, fed%Lz, fed%xin, fed%xiq, fed%wtq3, &
        T_au, &
        Venq, psi**2, fed%phihq, fed%Ap, fed%Aj, fed%Ax, fed%matd, &
        fed%spectral, &
        fed%phi_v, fed%jac_det, &
        Eh, Een, Ts, Exc, free_energy_, Hpsi=Hpsi)
!    print *, "Iteration:", iter
    psi_norm = integral(fed%nodes, fed%elems, fed%wtq3, psi**2)
!    print *, "Norm of psi:", psi_norm
!    print *, "mu =", mu
!    print *, "|ksi| =", sqrt(gamma_n)
!    print *, "theta =", theta
!    print *, "Summary of energies [a.u.]:"
!    print "('    Ts   = ', f14.8)", Ts
!    print "('    Een  = ', f14.8)", Een
!    print "('    Eee  = ', f14.8)", Eh
!    print "('    Exc  = ', f14.8)", Exc
!    print *, "   ---------------------"
    print "('    Etot = ', f14.8, ' a.u.')", free_energy_
    free_energies(iter) = free_energy_
    if (iter > 3) then
        last3 = maxval(free_energies(iter-3:iter)) - &
            minval(free_energies(iter-3:iter))
        if (last3 < energy_eps) then
            nq_pos = psi**2
            return
        end if
    end if
    Hpsi = Hpsi * 2*psi ! d/dpsi = 2 psi d/dn
    mu = 1._dp / Nelec * integral(fed%nodes, fed%elems, fed%wtq3, 0.5_dp * psi * Hpsi)
    ksi_prev = ksi
    ksi = 2*mu*psi - Hpsi
    gamma_n = max(integral(fed%nodes, fed%elems, fed%wtq3, ksi*(ksi-ksi_prev)), 0._dp)
    gamma_d = integral(fed%nodes, fed%elems, fed%wtq3, ksi_prev**2)
    phi = ksi + gamma_n / gamma_d * phi
    phi_prime = phi - 1._dp / Nelec *  integral(fed%nodes, fed%elems, fed%wtq3, phi * psi) * psi
    eta = sqrt(Nelec / integral(fed%nodes, fed%elems, fed%wtq3, phi_prime**2)) * phi_prime
end do
call stop_error("free_energy_minimization: The maximum number of iterations exceeded.")

contains

    real(dp) function func(theta) result(energy)
    real(dp), intent(in) :: theta
    psi_ = cos(theta) * psi + sin(theta) * eta
    call free_energy(fed%nodes, fed%elems, fed%in, fed%ib, fed%Nb, fed%Lx, fed%Ly, fed%Lz, fed%xin, fed%xiq, fed%wtq3, &
        T_au, &
        Venq, psi_**2, fed%phihq, fed%Ap, fed%Aj, fed%Ax, fed%matd, fed%spectral, &
        fed%phi_v, fed%jac_det, &
        Eh, Een, Ts, Exc, energy, Hpsi=Hpsi)
    end function

end subroutine



subroutine free_energy(nodes, elems, in, ib, Nb, Lx, Ly, Lz, xin, xiq, wtq3, &
    T_au, Venq, nq_pos, phihq, Ap, Aj, Ax, matd, spectral, &
    phi_v, jac_det, &
    Eh, Een, Ts, Exc, Etot, verbose, Hpsi)
real(dp), intent(in) :: nodes(:, :)
integer, intent(in) :: elems(:, :), in(:, :, :, :), ib(:, :, :, :), Nb
real(dp), intent(in) :: Lx, Ly, Lz, T_au
real(dp), intent(in) :: Venq(:, :, :, :), nq_pos(:, :, :, :), phihq(:, :), &
    xin(:), xiq(:), wtq3(:, :, :), &
    phi_v(:, :, :, :, :, :), jac_det, Ax(:)
integer, intent(in) :: Ap(:), Aj(:)
type(umfpack_numeric), intent(in) :: matd
real(dp), intent(out) :: Eh, Een, Ts, Exc, Etot
logical, intent(in) :: spectral
logical, intent(in), optional :: verbose
! If Hpsi is present, it
! returns "delta F / delta n", the functional derivative with respect to the
! density "n". Use the relation
!     d/dpsi = 2 psi d/dn
! to obtain the derivative with respect to psi (i.e. multiply Hpsi by 2*psi).
real(dp), intent(out) :: Hpsi(:, :, :, :)

real(dp), allocatable, dimension(:, :, :, :) :: y, F0, exc_density, &
    Vhq, nq_neutral, Vxc, dF0dn
real(dp), allocatable :: rhs(:), sol(:), fullsol(:)
real(dp) :: background, beta, tmp, dydn
integer :: p, i, j, k, m, Nn, Ne, Nq
logical :: verbose_
verbose_ = .false.
if (present(verbose)) verbose_ = verbose
Nn = size(nodes, 2)
Ne = size(elems, 2)
Nq = size(xiq)
p = size(xin) - 1
allocate(rhs(Nb), sol(Nb), fullsol(maxval(in)), Vhq(Nq, Nq, Nq, Ne))
allocate(y(Nq, Nq, Nq, Ne))
allocate(F0(Nq, Nq, Nq, Ne))
allocate(exc_density(Nq, Nq, Nq, Ne))
allocate(nq_neutral(Nq, Nq, Nq, Ne))
! Make the charge density net neutral (zero integral):
background = integral(nodes, elems, wtq3, nq_pos) / (Lx*Ly*Lz)
if (verbose_) then
    print *, "Total (positive) electronic charge: ", background * (Lx*Ly*Lz)
    print *, "Subtracting constant background (Q/V): ", background
end if
nq_neutral = nq_pos - background
if (verbose_) then
    print *, "Assembling RHS..."
end if
if (spectral) then
    call assemble_3d_coo_rhs_spectral(Ne, p, 4*pi*nq_neutral, jac_det, wtq3, &
        ib, rhs)
else
    call assemble_3d_coo_rhs(Ne, p, 4*pi*nq_neutral, jac_det, wtq3, &
        ib, phi_v, rhs)
end if
if (verbose_) then
    print *, "sum(rhs):    ", sum(rhs)
    print *, "integral rhs:", integral(nodes, elems, wtq3, nq_neutral)
    print *, "Solving..."
end if
if (WITH_UMFPACK) then
    call solve(Ap, Aj, Ax, sol, rhs, matd)
else
    sol = solve_cg(Ap, Aj, Ax, rhs, zeros(size(rhs)), 1e-12_dp, 800)
end if
if (verbose_) then
    print *, "Converting..."
end if
call c2fullc_3d(in, ib, sol, fullsol)
if (verbose_) then
    print *, "Transferring to quadrature points"
end if
if (spectral) then
    call fe2quad_3d_lobatto(elems, xiq, in, fullsol, Vhq)
else
    call fe2quad_3d(elems, xin, xiq, phihq, in, fullsol, Vhq)
end if
if (verbose_) then
    print *, "Done"
end if

! Hartree energy
Eh = integral(nodes, elems, wtq3, Vhq*nq_neutral) / 2
! Electron-nucleus energy
Een = integral(nodes, elems, wtq3, Venq*nq_neutral)
! Kinetic energy using Perrot parametrization
beta = 1/T_au
! The density must be positive, the f(y) fails for negative "y". Thus we use
! nq_pos.
y = pi**2 / sqrt(2._dp) * beta**(3._dp/2) * nq_pos
if (any(y < 0)) call stop_error("Density must be positive")
F0 = nq_pos / beta * f(y)
Ts = integral(nodes, elems, wtq3, F0)
! Exchange and correlation energy
do m = 1, Ne
do k = 1, Nq
do j = 1, Nq
do i = 1, Nq
    call xc_pz(nq_pos(i, j, k, m), exc_density(i, j, k, m), tmp)
end do
end do
end do
end do
Exc = integral(nodes, elems, wtq3, exc_density * nq_pos)
Etot = Ts + Een + Eh + Exc

! Calculate the derivative
allocate(Vxc(Nq, Nq, Nq, Ne), dF0dn(Nq, Nq, Nq, Ne))
dydn = pi**2 / sqrt(2._dp) * beta**(3._dp/2)
! F0 = nq_pos / beta * f(y)
! d F0 / d n =
dF0dn = 1 / beta * f(y) + nq_pos / beta * f(y, deriv=.true.) * dydn
! Exchange and correlation potential
do m = 1, Ne
do k = 1, Nq
do j = 1, Nq
do i = 1, Nq
    call xc_pz(nq_pos(i, j, k, m), tmp, Vxc(i, j, k, m))
end do
end do
end do
end do
Hpsi = dF0dn + Vhq + Venq + Vxc
Hpsi = Hpsi * Lx*Ly*Lz
end subroutine

subroutine free_energy2(Nelec, L, Nex, Ney, Nez, p, T_au, fnen, fn_pos, &
        Nq, quad_type, &
        Eh, Een, Ts, Exc, Nb)
! Higher level subroutine that does FE precalculation inside it, so it is slow,
! but the interface is simpler.
real(dp), intent(in) :: Nelec, L, T_au
integer, intent(in) :: p, Nq, quad_type, Nex, Ney, Nez
procedure(func_xyz) :: fnen ! (negative) ionic particle density
procedure(func_xyz) :: fn_pos ! (positive) electronic particle density
real(dp), intent(out) :: Eh, Een, Ts, Exc
integer, intent(out) :: Nb

real(dp), allocatable, dimension(:, :, :, :) :: nenq_pos, nq_pos
type(fe_data) :: fed

call initialize_fe(L, Nex, Ney, Nez, p, Nq, quad_type, fed)

allocate(nenq_pos(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(nq_pos(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
nenq_pos = func2quad(fed%nodes, fed%elems, fed%xiq, fnen)
nq_pos = func2quad(fed%nodes, fed%elems, fed%xiq, fn_pos)

call free_energy2_low_level(Nelec, T_au, nenq_pos, nq_pos, &
        fed, Eh, Een, Ts, Exc)

Nb = fed%Nb
end subroutine

subroutine free_energy2_low_level(Nelec, T_au, nenq_pos, nq_pos, &
        ! Geometry, finite element specification and internal datastructures
        fed, &
        ! Output arguments
        Eh, Een, Ts, Exc)
! Lower level subroutine that doesn't do any FE precalculation, so it is fast.
real(dp), intent(in) :: Nelec
real(dp), intent(in) :: nenq_pos(:, :, :, :)
real(dp), intent(inout) :: nq_pos(:, :, :, :)
real(dp), intent(in) :: T_au
type(fe_data), intent(in) :: fed
real(dp), intent(out) :: Eh, Een, Ts, Exc

real(dp), allocatable, dimension(:, :, :, :) :: Hpsi, &
    psi, Venq, &
    nenq_neutral
integer :: max_iter
real(dp) :: brent_eps, free_energy_
real(dp) :: psi_norm
real(dp) :: background
real(dp), allocatable :: rhs(:), sol(:), fullsol(:), fullsol2(:)

brent_eps = 1e-3_dp
max_iter = 200

allocate(nenq_neutral(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(Venq(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(Hpsi(fed%Nq, fed%Nq, fed%Nq, fed%Ne))
allocate(psi(fed%Nq, fed%Nq, fed%Nq, fed%Ne))


! Calculate Venq
allocate(rhs(fed%Nb), sol(fed%Nb), fullsol(maxval(fed%in)))
allocate(fullsol2(maxval(fed%in)))
background = integral(fed%nodes, fed%elems, fed%wtq3, nenq_pos) / &
    (fed%Lx*fed%Ly*fed%Lz)
print *, "Total (negative) ionic charge: ", background * (fed%Lx*fed%Ly*fed%Lz)
print *, "Subtracting constant background (Q/V): ", background
nenq_neutral = nenq_pos - background
print *, "Assembling RHS..."
if (fed%spectral) then
    call assemble_3d_coo_rhs_spectral(fed%Ne, fed%p, 4*pi*nenq_neutral, &
        fed%jac_det, fed%wtq3, fed%ib, rhs)
else
    call assemble_3d_coo_rhs(fed%Ne, fed%p, 4*pi*nenq_neutral, &
        fed%jac_det, fed%wtq3, fed%ib, fed%phi_v, rhs)
end if
print *, "sum(rhs):    ", sum(rhs)
print *, "integral rhs:", integral(fed%nodes, fed%elems, fed%wtq3, nenq_neutral)
print *, "Solving..."
if (WITH_UMFPACK) then
    call solve(fed%Ap, fed%Aj, fed%Ax, sol, rhs, fed%matd)
else
    sol = solve_cg(fed%Ap, fed%Aj, fed%Ax, rhs, zeros(size(rhs)), 1e-12_dp, 800)
end if
print *, "Converting..."
call c2fullc_3d(fed%in, fed%ib, sol, fullsol)
if (fed%spectral) then
    call fe2quad_3d_lobatto(fed%elems, fed%xiq, fed%in, fullsol, Venq)
else
    call fe2quad_3d(fed%elems, fed%xin, fed%xiq, fed%phihq, fed%in, fullsol, Venq)
end if
print *, "Done"


psi = sqrt(nq_pos)
psi_norm = integral(fed%nodes, fed%elems, fed%wtq3, psi**2)
print *, "Initial norm of psi:", psi_norm
psi = sqrt(Nelec / psi_norm) * psi
psi_norm = integral(fed%nodes, fed%elems, fed%wtq3, psi**2)
print *, "norm of psi:", psi_norm
! This returns H[n] = delta F / delta n, we save it to the Hpsi variable to
! save space:
call free_energy(fed%nodes, fed%elems, fed%in, fed%ib, fed%Nb, fed%Lx, fed%Ly, fed%Lz, fed%xin, fed%xiq, fed%wtq3, T_au, &
    Venq, psi**2, fed%phihq, fed%Ap, fed%Aj, fed%Ax, fed%matd, fed%spectral, &
    fed%phi_v, fed%jac_det, &
    Eh, Een, Ts, Exc, free_energy_, Hpsi=Hpsi)
end subroutine


subroutine radial_density_fourier(R, V, L, Z, Ng, Rnew, density)
real(dp), intent(in) :: R(:), V(:), L, Z, Rnew(:)
real(dp), intent(out) :: density(:)
real(dp), allocatable :: density_ft(:), w(:)
real(dp) :: Rp(size(R))
integer, intent(in) :: Ng ! Number of PW
integer :: j
real(dp) :: dk, Rc
! Rp is the derivative of the mesh R'(t), which for uniform mesh is equal to
! the mesh step (rmax-rmin)/N:
Rp = (R(size(R)) - R(1)) / (size(R)-1)
Rc = R(size(R))
dk = 2*pi/L
! To fill out a 3D grid, we would use this:
!allocate(Vk(3*(Ng/2+1)**2))
!do j = 1, 3*(Ng/2+1)**2
!    w = sqrt(real(j, dp))*dk
!    Vk(j) = 4*pi/(L**3 * w**2) * (w*integrate(Rp, R*sin(w*R)*V) + cos(w*Rc))
!end do

! We can use denser mesh, but this looks sufficient:
allocate(density_ft(Ng), w(Ng))
forall(j=1:Ng) w(j) = (j-1)*dk
forall (j=1:Ng)
    density_ft(j) = w(j)*integrate(Rp, R*sin(w(j)*R)*V) + cos(w(j)*Rc)
end forall
density_ft = Z*density_ft
! The potential (infinite for w=0) would be calculated using:
!Vk = 4*pi/w**2 * density_ft

forall (j=1:size(Rnew))
    density(j) = sum(dk * w**2*sinc(w*Rnew(j))*density_ft) / (2*pi**2)
end forall
end subroutine

real(dp) pure function integrate(Rp, f) result(s)
real(dp), intent(in) :: Rp(:), f(:)
! Choose one from the integration rules below:
s = integrate_trapz_1(Rp, f)
!s = integrate_trapz_3(Rp, f)
!s = integrate_trapz_5(Rp, f)
!s = integrate_trapz_7(Rp, f)
!s = integrate_simpson(Rp, f)
!s = integrate_adams(Rp, f)
end function

real(dp) elemental function sinc(x) result(r)
real(dp), intent(in) :: x
if (abs(x) < 1e-8_dp) then
    r = 1
else
    r = sin(x)/x
end if
end function


end module
