module ofmd_utils

use types, only: dp
use feutils, only: phih, dphih
use fe_mesh, only: cartesian_mesh_3d, define_connect_tensor_3d, &
    c2fullc_3d, fe2quad_3d, vtk_save, fe_eval_xyz, line_save
use poisson3d_assembly, only: assemble_3d, integral, func2quad, func_xyz
use feutils, only: get_parent_nodes, get_parent_quad_pts_wts
use linalg, only: solve
use isolve, only: solve_cg
use utils, only: assert, zeros, stop_error
use constants, only: pi
use xc, only: xc_pz
use optimize, only: bracket, brent
implicit none
private
public free_energy_min, read_pseudo

contains

subroutine free_energy_min(L, Nex, Ney, Nez, p, T_au, fnen, fn_pos, Eh, Een, &
    Ts, Exc, Nb)
integer, intent(in) :: p
procedure(func_xyz) :: fnen ! (negative) ionic particle density
procedure(func_xyz) :: fn_pos ! (positive) electronic particle density
real(dp), intent(in) :: L, T_au
real(dp), intent(out) :: Eh, Een, Ts, Exc
integer, intent(out) :: Nb

integer :: Nn, Ne, ibc
! nodes(:, i) are the (x,y) coordinates of the i-th mesh node
real(dp), allocatable :: nodes(:, :)
integer, allocatable :: elems(:, :) ! elems(:, i) are nodes of the i-th element
integer :: Nq
real(dp), allocatable :: xin(:), xiq(:), wtq(:), &
        wtq3(:, :, :), phihq(:, :), dphihq(:, :), free_energies(:)
real(dp), allocatable, dimension(:, :, :, :) :: nenq_pos, nq_pos, Hpsi, &
    psi, psi_, psi_prev, ksi, ksi_prev, phi, phi_prime, eta
integer, allocatable :: in(:, :, :, :), ib(:, :, :, :)
integer :: i, j, k, iter, max_iter
integer, intent(in) :: Nex, Ney, Nez
real(dp) :: Lx, Ly, Lz, mu, energy_eps, last3, brent_eps, free_energy_, &
    gamma_d, gamma_n, theta, theta_a, theta_b, theta_c

energy_eps = 1e-10_dp
brent_eps = 1e-10_dp
max_iter = 30

ibc = 3 ! Periodic boundary condition

Lx = L
Ly = L
Lz = L

call cartesian_mesh_3d(Nex, Ney, Nez, &
    [-Lx/2, -Ly/2, -Lz/2], [Lx/2, Ly/2, Lz/2], nodes, elems)
Nn = size(nodes, 2)
Ne = size(elems, 2)
Nq = 20

print *, "Number of nodes:", Nn
print *, "Number of elements:", Ne
print *, "Nq =", Nq
print *, "p =", p
allocate(xin(p+1))
call get_parent_nodes(2, p, xin)
allocate(xiq(Nq), wtq(Nq), wtq3(Nq, Nq, Nq))
call get_parent_quad_pts_wts(1, Nq, xiq, wtq)
forall(i=1:Nq, j=1:Nq, k=1:Nq) wtq3(i, j, k) = wtq(i)*wtq(j)*wtq(k)
allocate(phihq(size(xiq), size(xin)))
allocate(dphihq(size(xiq), size(xin)))
! Tabulate parent basis at quadrature points
forall(i=1:size(xiq), j=1:size(xin))  phihq(i, j) =  phih(xin, j, xiq(i))
forall(i=1:size(xiq), j=1:size(xin)) dphihq(i, j) = dphih(xin, j, xiq(i))

call define_connect_tensor_3d(Nex, Ney, Nez, p, 1, in)
call define_connect_tensor_3d(Nex, Ney, Nez, p, ibc, ib)
Nb = maxval(ib)
print *, "DOFs =", Nb
allocate(nenq_pos(Nq, Nq, Nq, Ne))
allocate(nq_pos(Nq, Nq, Nq, Ne))
allocate(Hpsi(Nq, Nq, Nq, Ne))
allocate(psi(Nq, Nq, Nq, Ne))
allocate(psi_(Nq, Nq, Nq, Ne))
allocate(psi_prev(Nq, Nq, Nq, Ne))
allocate(phi(Nq, Nq, Nq, Ne))
allocate(phi_prime(Nq, Nq, Nq, Ne))
allocate(ksi(Nq, Nq, Nq, Ne))
allocate(ksi_prev(Nq, Nq, Nq, Ne))
allocate(eta(Nq, Nq, Nq, Ne))

nenq_pos = func2quad(nodes, elems, xiq, fnen)
nq_pos = func2quad(nodes, elems, xiq, fn_pos)

psi = sqrt(nq_pos)
call free_energy_derivative(nodes, elems, in, ib, Nb, Lx, Ly, Lz, xin, &
    xiq, wtq3, T_au, &
    nenq_pos, psi**2, phihq, dphihq, Hpsi) ! Hpsi = <x|H|psi> = H psi(x, y, z)
Hpsi = Hpsi * 2*psi ! d/dpsi = 2 psi d/dn
mu = 1._dp / Ne * integral(nodes, elems, wtq3, 0.5_dp * psi * Hpsi)
ksi = 2*mu*psi - Hpsi
phi = ksi
phi_prime = phi - 1._dp / Ne *  integral(nodes, elems, wtq3, phi * psi) * psi
eta = sqrt(Ne / integral(nodes, elems, wtq3, phi_prime**2)) * phi_prime
theta = pi/2
call free_energy(nodes, elems, in, ib, Nb, Lx, Ly, Lz, xin, xiq, wtq3, T_au, &
    nenq_pos, psi**2, phihq, dphihq, Eh, Een, Ts, Exc, free_energy_)
allocate(free_energies(max_iter))
do iter = 1, max_iter
    theta_a = 0
    theta_b = mod(theta, 2*pi)
    call bracket(f, theta_a, theta_b, theta_c, 100._dp, 20, verbose=.true.)
    print *, "bracket:", theta_a, theta_b, theta_c
    call brent(f, theta_a, theta_b, theta_c, brent_eps, 50, theta, &
        free_energy_, verbose=.true.)
    ! TODO: We probably don't need to recalculate free_energy_ here:
    psi_prev = psi
    psi = cos(theta) * psi + sin(theta) * eta
    call free_energy(nodes, elems, in, ib, Nb, Lx, Ly, Lz, xin, xiq, wtq3, &
        T_au, &
        nenq_pos, psi**2, phihq, dphihq, Eh, Een, Ts, Exc, free_energy_)
    print *, "Summary of energies [a.u.]:"
    print "('    Ts   = ', f14.8)", Ts
    print "('    Een  = ', f14.8)", Een
    print "('    Eee  = ', f14.8)", Eh
    print "('    Exc  = ', f14.8)", Exc
    print *, "   ---------------------"
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
    ! TODO: calculate the derivative together with the free energy, as pretty
    ! much the whole calculation is shared.
    call free_energy_derivative(nodes, elems, in, ib, Nb, Lx, Ly, Lz, xin, &
        xiq, wtq3, T_au, &
        nenq_pos, psi**2, phihq, dphihq, Hpsi)
    Hpsi = Hpsi * 2*psi ! d/dpsi = 2 psi d/dn
    ksi_prev = ksi
    ksi = 2*mu*psi - Hpsi
    gamma_n = integral(nodes, elems, wtq3, ksi**2)
    gamma_d = integral(nodes, elems, wtq3, ksi_prev**2)
    phi = ksi + gamma_n / gamma_d * phi
    phi_prime = phi - 1._dp / Ne *  integral(nodes, elems, wtq3, phi * psi) * psi
    eta = sqrt(Ne / integral(nodes, elems, wtq3, phi_prime**2)) * phi_prime
end do
call stop_error("free_energy_minimization: The maximum number of iterations exceeded.")

contains

    real(dp) function f(theta) result(energy)
    real(dp), intent(in) :: theta
    psi_ = cos(theta) * psi + sin(theta) * eta
    call free_energy(nodes, elems, in, ib, Nb, Lx, Ly, Lz, xin, xiq, wtq3, &
        T_au, &
        nenq_pos, psi_**2, phihq, dphihq, Eh, Een, Ts, Exc, energy)
    end function

end subroutine



subroutine free_energy(nodes, elems, in, ib, Nb, Lx, Ly, Lz, xin, xiq, wtq3, &
    T_au, nenq_pos, nq_pos, phihq, dphihq, Eh, Een, Ts, Exc, Etot)
real(dp), intent(in) :: nodes(:, :)
integer, intent(in) :: elems(:, :), in(:, :, :, :), ib(:, :, :, :), Nb
real(dp), intent(in) :: Lx, Ly, Lz, T_au
real(dp), intent(in) :: nenq_pos(:, :, :, :), nq_pos(:, :, :, :), phihq(:, :), &
    dphihq(:, :), xin(:), xiq(:), wtq3(:, :, :)
real(dp), intent(out) :: Eh, Een, Ts, Exc, Etot

real(dp), allocatable, dimension(:, :, :, :) :: y, F0, exc_density, &
    nq_neutral, Venq, Vhq, nenq_neutral
integer, allocatable :: Ap(:), Aj(:)
real(dp), allocatable :: Ax(:), rhs(:), sol(:), fullsol(:)
real(dp) :: background, beta, tmp
integer :: i, j, k, m, Nn, Ne, Nq
Nn = size(nodes, 2)
Ne = size(elems, 2)
Nq = size(xiq)
allocate(rhs(Nb), sol(Nb), fullsol(maxval(in)), Vhq(Nq, Nq, Nq, Ne))
allocate(y(Nq, Nq, Nq, Ne))
allocate(F0(Nq, Nq, Nq, Ne))
allocate(exc_density(Nq, Nq, Nq, Ne))
allocate(nq_neutral(Nq, Nq, Nq, Ne))
allocate(Venq(Nq, Nq, Nq, Ne))
! Make the charge density net neutral (zero integral):
background = integral(nodes, elems, wtq3, nq_pos) / (Lx*Ly*Lz)
print *, "Total (positive) electronic charge: ", background * (Lx*Ly*Lz)
print *, "Subtracting constant background (Q/V): ", background
nq_neutral = nq_pos - background
call assemble_3d(xin, nodes, elems, ib, xiq, wtq3, phihq, dphihq, &
    4*pi*nq_neutral, Ap, Aj, Ax, rhs)
print *, "sum(rhs):    ", sum(rhs)
print *, "integral rhs:", integral(nodes, elems, wtq3, nq_neutral)
print *, "Solving..."
!sol = solve(A, rhs)
sol = solve_cg(Ap, Aj, Ax, rhs, zeros(size(rhs)), 1e-12_dp, 400)
call c2fullc_3d(in, ib, sol, fullsol)
call fe2quad_3d(elems, xin, xiq, phihq, in, fullsol, Vhq)

background = integral(nodes, elems, wtq3, nenq_pos) / (Lx*Ly*Lz)
print *, "Total (negative) ionic charge: ", background * (Lx*Ly*Lz)
print *, "Subtracting constant background (Q/V): ", background
nenq_neutral = nenq_pos - background
call assemble_3d(xin, nodes, elems, ib, xiq, wtq3, phihq, dphihq, &
    4*pi*nenq_neutral, Ap, Aj, Ax, rhs)
print *, "sum(rhs):    ", sum(rhs)
print *, "integral rhs:", integral(nodes, elems, wtq3, nenq_neutral)
print *, "Solving..."
sol = solve_cg(Ap, Aj, Ax, rhs, zeros(size(rhs)), 1e-12_dp, 400)
call c2fullc_3d(in, ib, sol, fullsol)
call fe2quad_3d(elems, xin, xiq, phihq, in, fullsol, Venq)

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
end subroutine

subroutine free_energy_derivative(nodes, elems, in, ib, Nb, Lx, Ly, Lz, xin, &
        xiq, wtq3, T_au, &
        nenq_pos, nq_pos, phihq, dphihq, Hpsi)
! Returns "delta F / delta n", the functional derivative with respect to the
! density "n". Use the relation
!     d/dpsi = 2 psi d/dn
! to obtain the derivative with respect to psi (i.e. multiply Hpsi by 2*psi).
real(dp), intent(in) :: nodes(:, :)
integer, intent(in) :: elems(:, :), in(:, :, :, :), ib(:, :, :, :), Nb
real(dp), intent(in) :: Lx, Ly, Lz, T_au
real(dp), intent(in) :: nenq_pos(:, :, :, :), nq_pos(:, :, :, :), phihq(:, :), &
    dphihq(:, :), xin(:), xiq(:), wtq3(:, :, :)
real(dp), intent(out) :: Hpsi(:, :, :, :)

real(dp), allocatable, dimension(:, :, :, :) :: y, dF0dn, exc_density, &
    nq_neutral, Venq, Vhq, nenq_neutral, Vxc
integer, allocatable :: Ap(:), Aj(:)
real(dp), allocatable :: Ax(:), rhs(:), sol(:), fullsol(:)
real(dp) :: background, beta, tmp, dydn
integer :: Nn, Ne, Nq, i, j, k, m
Nn = size(nodes, 2)
Ne = size(elems, 2)
Nq = size(xiq)
allocate(rhs(Nb), sol(Nb), fullsol(maxval(in)), Vhq(Nq, Nq, Nq, Ne))
allocate(y(Nq, Nq, Nq, Ne))
allocate(dF0dn(Nq, Nq, Nq, Ne))
allocate(exc_density(Nq, Nq, Nq, Ne))
allocate(nq_neutral(Nq, Nq, Nq, Ne))
allocate(Venq(Nq, Nq, Nq, Ne))
allocate(Vxc(Nq, Nq, Nq, Ne))
! Make the charge density net neutral (zero integral):
background = integral(nodes, elems, wtq3, nq_pos) / (Lx*Ly*Lz)
print *, "Total (positive) electronic charge: ", background * (Lx*Ly*Lz)
print *, "Subtracting constant background (Q/V): ", background
nq_neutral = nq_pos - background
call assemble_3d(xin, nodes, elems, ib, xiq, wtq3, phihq, dphihq, &
    4*pi*nq_neutral, Ap, Aj, Ax, rhs)
print *, "sum(rhs):    ", sum(rhs)
print *, "integral rhs:", integral(nodes, elems, wtq3, nq_neutral)
print *, "Solving..."
!sol = solve(A, rhs)
sol = solve_cg(Ap, Aj, Ax, rhs, zeros(size(rhs)), 1e-12_dp, 400)
call c2fullc_3d(in, ib, sol, fullsol)
call fe2quad_3d(elems, xin, xiq, phihq, in, fullsol, Vhq)

background = integral(nodes, elems, wtq3, nenq_pos) / (Lx*Ly*Lz)
print *, "Total (negative) ionic charge: ", background * (Lx*Ly*Lz)
print *, "Subtracting constant background (Q/V): ", background
nenq_neutral = nenq_pos - background
call assemble_3d(xin, nodes, elems, ib, xiq, wtq3, phihq, dphihq, &
    4*pi*nenq_neutral, Ap, Aj, Ax, rhs)
print *, "sum(rhs):    ", sum(rhs)
print *, "integral rhs:", integral(nodes, elems, wtq3, nenq_neutral)
print *, "Solving..."
sol = solve_cg(Ap, Aj, Ax, rhs, zeros(size(rhs)), 1e-12_dp, 400)
call c2fullc_3d(in, ib, sol, fullsol)
call fe2quad_3d(elems, xin, xiq, phihq, in, fullsol, Venq)

beta = 1/T_au
y    = pi**2 / sqrt(2._dp) * beta**(3._dp/2) * nq_pos
dydn = pi**2 / sqrt(2._dp) * beta**(3._dp/2)
if (any(y < 0)) call stop_error("Density must be positive")
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
end subroutine

subroutine read_pseudo(filename, R, V, Z, Ediff)
! Reads the pseudopotential from the file 'filename'.
character(len=*), intent(in) :: filename   ! File to read from, e.g. "H.pseudo"
real(dp), allocatable, intent(out) :: R(:) ! radial grid [0, Rcut]
! potential on the radial grid. The potential smoothly changes into -1/R for
! r > Rcut, where Rcut = R(size(R)) is the cut-off radius
real(dp), allocatable, intent(out) :: V(:)
real(dp), intent(out) :: Z     ! Nuclear charge
real(dp), intent(out) :: Ediff ! The energy correction
real(dp) :: Rcut
integer :: N, i, u
open(newunit=u, file=filename, status="old")
read(u, *) Z, N, Rcut, Ediff
allocate(R(N-1), V(N-1))
! The first potential value is zero in the file, so we skip it
read(u, *) R(1), V(1)
do i = 1, N-1
    read(u, *) R(i), V(i)
end do
close(u)
! The file contains a grid from [0, 1], so we need to rescale it:
R = R*Rcut
! We need to add the minus sign to the potential ourselves:
V = -V
end subroutine

real(dp) elemental function f(y, deriv)
! Function f(y) from Appendix A in [1].
!
! [1] Perrot, F. (1979). Gradient correction to the statistical electronic free
! energy at nonzero temperatures: Application to equation-of-state
! calculations. Physical Review A, 20(2), 586–594.
real(dp), intent(in) :: y ! must be positive
! if deriv == .true. compute df/dy instead. Default .false.
logical, intent(in), optional :: deriv
real(dp), parameter :: y0 = 3*pi/(4*sqrt(2._dp))
real(dp), parameter :: c(*) = [-0.8791880215_dp, 0.1989718742_dp, &
    0.1068697043e-2_dp, -0.8812685726e-2_dp, 0.1272183027e-1_dp, &
    -0.9772758583e-2_dp, 0.3820630477e-2_dp, -0.5971217041e-3_dp]
real(dp), parameter :: d(*) = [0.7862224183_dp, -0.1882979454e1_dp, &
    0.5321952681_dp, 0.2304457955e1_dp, -0.1614280772e2_dp, &
    0.5228431386e2_dp, -0.9592645619e2_dp, 0.9462230172e2_dp, &
    -0.3893753937e2_dp]
real(dp) :: u
integer :: i
logical :: deriv_
deriv_ = .false.
if (present(deriv)) deriv_ = deriv

if (.not. deriv_) then
    if (y <= y0) then
        f = log(y)
        do i = 0, 7
            f = f + c(i+1) * y**i
        end do
    else
        u = y**(2._dp / 3)
        f = 0
        do i = 0, 8
            f = f + d(i+1) / u**(2*i-1)
        end do
        ! Note: Few terms in [1] have "y" instead of "u" in them for y > y0, but
        ! that is obviously a typo.
    end if
else
    if (y <= y0) then
        f = 1 / y
        do i = 0, 6
            f = f + (i+1) * c(i+2) * y**i
        end do
    else
        u = y**(2._dp / 3)
        f = 0
        do i = 0, 8
            f = f - (2*i-1) * d(i+1) / u**(2*i)
        end do
        f = f * 2._dp/3 / y**(1._dp/3)
    end if
end if
end function

end module


! ------------------------------------------------------------------------


program ofmd
use types, only: dp
use ofmd_utils, only: free_energy_min, read_pseudo
use constants, only: Ha2eV, pi
use utils, only: loadtxt
use splines, only: spline3pars, iixmin, poly3
use interp3d, only: trilinear
implicit none
real(dp) :: Eh, Een, Ts, Exc, Etot
integer :: p, DOF !, Nx, Ny, Nz, u
real(dp) :: Z, Ediff
real(dp), allocatable :: R(:), V(:), D(:, :), c(:, :), values(:, :, :)
real(dp) :: Rcut, L, T_eV, T_au

!open(newunit=u, file="plots/Ven_reg128.txt", status="old")
!read(u, *) Nx, Ny, Nz
!allocate(values(Nx, Ny, Nz))
!read(u, *) values
!close(u)
!call read_pseudo("H.pseudo", R, V, Z, Ediff)
allocate(R(1), V(1)); Z=1; Ediff=0
!call loadtxt("Venr.txt", D)
!allocate(c(0:4, size(D, 1)-1))
!call spline3pars(D(:, 1), D(:, 2), [2, 2], [0._dp, 0._dp], c)
Rcut = R(size(R))
Rcut = 0.3_dp
p = 5
L = 2
T_eV = 0.0862_dp
T_au = T_ev / Ha2eV
call free_energy_min(L, 5, 5, 5, p, T_au, nen, ne, Eh, Een, Ts, Exc, DOF)
Etot = Ts + Een + Eh + Exc
print *, "p =", p
print *, "DOF =", DOF
print *, "Rcut =", Rcut
print *, "T_au =", T_au
print *, "Summary of energies [a.u.]:"
print "('    Ts   = ', f14.8)", Ts
print "('    Een  = ', f14.8)", Een
print "('    Eee  = ', f14.8)", Eh
print "('    Exc  = ', f14.8)", Exc
print *, "   ---------------------"
print "('    Etot = ', f14.8, ' a.u. = ', f14.8, ' eV')", Etot, Etot*Ha2eV

contains

real(dp) function nen(x, y, z_) result(n)
real(dp), intent(in) :: x, y, z_
real(dp), parameter :: alpha = 12
real(dp) :: r
r = sqrt(x**2+y**2+z_**2)
! This density:
n = -Z*alpha**3/pi**(3._dp/2)*exp(-alpha**2*R**2)
! Corresponds to the potential:
!V = -Z*erf(alpha*R)/R
end function

real(dp) function Ven_splines(x_, y_, z_) result(V)
real(dp), intent(in) :: x_, y_, z_
real(dp) :: r
integer :: ip
! One atom in the center:
r = sqrt(x_**2+y_**2+z_**2)
if (r >= 1) r = 1
ip = 0
ip = iixmin(r, D(:, 1), ip)
V = poly3(r, c(:, ip))
end function

real(dp) function Ven_interp3d(x, y, z) result(V)
real(dp), intent(in) :: x, y, z
V = trilinear([x, y, z], [-L/2, -L/2, -L/2], [L/2, L/2, L/2], &
        values)
end function

real(dp) function ne(x, y, z) result(n)
real(dp), intent(in) :: x, y, z
real(dp), parameter :: alpha = 5, Z_ = 1
real(dp) :: r
r = sqrt(x**2+y**2+z**2)
n = Z_*alpha**3/pi**(3._dp/2)*exp(-alpha**2*R**2)
end function

end program
