! Configure, compile and run with:
! CC=mpicc CXX=mpic++ FC=mpif90 cmake -DCOMMON_DIR=$HASHSTACK -DWITH_BDDC=yes .
! make
! cd src/tests/fem
! mpirun -np 2 ./test_free_energy_bddc
!
! Conejo:
! module load gcc/4.8.2 openmpi
! export HASHSTACK=/users/certik/repos/hashstack/default
! CC=mpicc CXX=mpic++ FC=mpif90 cmake -DCOMMON_DIR=$HASHSTACK -DWITH_BDDC=yes .
! make
! cd src/tests/fem
! mpiexec -n 2 ./test_free_energy_bddc

module test_free_energy_bddc_utils

use types, only: dp
use feutils, only: phih, dphih
use fe_mesh, only: cartesian_mesh_3d, define_connect_tensor_3d, &
    c2fullc_3d, fe2quad_3d_lobatto, vtk_save, fe_eval_xyz, line_save, &
    cartesian_mesh_3d_mask
use poisson3d_assembly, only: assemble_3d, integral, func2quad, func_xyz, &
    assemble_3d_precalc, assemble_3d_coo_A_subdomain_spectral, &
    assemble_3d_coo_rhs_spectral
use feutils, only: get_parent_nodes, get_parent_quad_pts_wts
use linalg, only: solve
use isolve, only: solve_cg
use utils, only: assert, zeros, stop_error
use constants, only: pi
use xc, only: xc_pz
use mpi_interface, only: mpi_comm_rank, mpi_double_precision, mpi_comm_size, &
    mpi_bcast
use bddcml_interface_fortran, only: bddcml_init, &
    bddcml_upload_subdomain_data, bddcml_setup_preconditioner, bddcml_solve, &
    bddcml_download_global_solution, bddcml_change_subdomain_data, &
    bddcml_setup_new_data, bddcml_finalize
implicit none
private
public free_energy, read_pseudo

contains

subroutine free_energy(comm, L, Nex, Ney, Nez, nsubx, nsuby, nsubz, p, T_au, fnen, fn_pos, Eh, Een, Ts, &
        Exc, Nb)
integer, intent(in) :: comm
integer, intent(in) :: p
procedure(func_xyz) :: fnen ! (negative) ionic particle density
procedure(func_xyz) :: fn_pos ! (positive) electronic particle density
real(dp), intent(in) :: L, T_au
real(dp), intent(out) :: Eh, Een, Ts, Exc
integer, intent(out) :: Nb

integer :: Nn, Ne, ibc
! nodes(:, i) are the (x,y,z) coordinates of the i-th mesh node
real(dp), allocatable :: nodes(:, :)
integer, allocatable :: elems(:, :) ! elems(:, i) are nodes of the i-th element
integer, allocatable :: mask_elems(:)
integer :: Nq
real(dp), allocatable :: xin(:), xiq(:), wtq(:), &
        rhs(:), sol(:), &
        fullsol(:), Vhq(:, :, :, :), wtq3(:, :, :), phihq(:, :), dphihq(:, :),&
        nenq(:, :, :, :), &
        Venq(:, :, :, :), y(:, :, :, :), F0(:, :, :, :), &
        exc_density(:, :, :, :), nq_pos(:, :, :, :), nq_neutral(:, :, :, :)
integer, allocatable :: in(:, :, :, :), ib(:, :, :, :) !Ap(:), Aj(:)
integer :: i, j, k, m
integer, intent(in) :: Nex, Ney, Nez
integer, intent(in) :: nsubx, nsuby, nsubz
real(dp) :: background
real(dp) :: Lx, Ly, Lz
real(dp) :: beta, tmp
integer :: myid, ierr, nproc
!real(dp) :: lx, ly, lz
real(dp) :: jac_det
integer, allocatable :: matAi(:), matAj(:)
real(dp),  allocatable :: matAx(:)
integer :: idx
integer :: Nesub ! number of elements on a subdomain
integer :: Nbelem
integer, allocatable :: global_to_local(:), local_to_global(:)
integer :: Nbsub ! number of basis functions on a subdomain
integer, allocatable :: elems_sub(:), isegns(:), nnets(:), nndfs(:)
integer, allocatable :: ifixs(:)
real(dp), allocatable :: xyzs(:, :), xyz_global(:, :), fixvs(:), rhss(:), &
    sols(:)
real(dp) :: user_constraints(1), element_data(1), dof_data(1)
real(dp), dimension(p+1) :: xp, yp, zp
real(dp) :: elx, ely, elz
integer :: iax, iay, iaz
integer :: e
real(dp) :: jacx, jacy, jacz
integer :: nlevels, nsub_loc_1, verbose_level
integer :: nsublev(2)
integer :: is_rhs_complete_int
integer :: is_assembled_int
integer :: matrixtype

real(dp) :: condition_number
integer :: num_iter, converged_reason
integer :: Asize

! *************************
! KRYLOV METHOD PARAMETERS:
! *************************

! Krylov subspace iterative method to be used
!     -1 - use solver defaults
!     0 - PCG
!     1 - BICGSTAB (choose for general symmetric and general matrices)
      integer,parameter :: krylov_method = 0

! use recycling of Krylov subspace
!     0 - no recycling used
!     1 - basis of the Krylov subspace will be orthogonalized and also used for new right hand sides
      integer,parameter :: recycling_int = 0
! size of the Krylov subspace basis to store
      integer,parameter :: max_number_of_stored_vectors = 50

! maximum number of iterations of a Krylov subspace method
      integer,parameter :: maxit = 500

! maximum number of iterations of a Krylov subspace method with non-decreasing residual
      integer,parameter :: ndecrmax = 50

! relative precision of the Krylov subspace method ||residual||/||right-hand side||
      real(dp),parameter :: tol = 1e-12_dp

! *******************************
! BDDC PRECONDITIONER PARAMETERS:
! *******************************


integer,parameter :: use_preconditioner_defaults = 0

! use arithmetic constraints on edges and faces?
integer,parameter :: use_arithmetic_constraints = 1

! use adaptive constraints on faces?
integer,parameter :: use_adaptive_constraints = 0

! use user constraints? - not used in this example
integer,parameter :: use_user_constraints = 0

! what type of weights use on interface?
! 0 - weights by cardinality
! 1 - weights by diagonal stiffness
! 2 - weights based on first row of element data
! 3 - weights based on dof data
! 4 - weights by Marta Certikova - unit load
! 5 - weights by Marta Certikova - unit jump
! 6 - weights by Schur row sums for whole subdomain
! 7 - weights by Schur row sums computed face by face
integer,parameter :: weights_type = 0

! should parallel division be used (ParMETIS instead of METIS) on the first level?
integer,parameter :: parallel_division = 1

integer :: isub

call MPI_COMM_RANK(comm,myid,ierr)
call MPI_COMM_SIZE(comm,nproc,ierr)
if (nsubx*nsuby*nsubz /= nproc) call stop_error("nproc must be equal to the number of subdomains")

ibc = 3 ! Periodic boundary condition

Lx = L
Ly = L
Lz = L

call cartesian_mesh_3d(Nex, Ney, Nez, &
    [-Lx/2, -Ly/2, -Lz/2], [Lx/2, Ly/2, Lz/2], nodes, elems)
call cartesian_mesh_3d_mask(myid, Nex, Ney, Nez, nsubx, nsuby, nsubz, &
        mask_elems)
!print *, "myid = ", myid, "mask_elems = ", mask_elems
!call MPI_BARRIER(comm, ierr)
!stop "OK"
Nn = size(nodes, 2)
Ne = size(elems, 2)
Nq = p+1

print *, "Number of nodes:", Nn
print *, "Number of elements:", Ne
print *, "Nq =", Nq
print *, "p =", p
allocate(xin(p+1))
call get_parent_nodes(2, p, xin)
allocate(xiq(Nq), wtq(Nq), wtq3(Nq, Nq, Nq))
call get_parent_quad_pts_wts(2, Nq, xiq, wtq)
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
allocate(rhs(Nb), sol(Nb), fullsol(maxval(in)), Vhq(Nq, Nq, Nq, Ne))
allocate(Venq(Nq, Nq, Nq, Ne))
allocate(nenq(Nq, Nq, Nq, Ne))
allocate(y(Nq, Nq, Nq, Ne))
allocate(F0(Nq, Nq, Nq, Ne))
allocate(exc_density(Nq, Nq, Nq, Ne))
allocate(nq_pos(Nq, Nq, Nq, Ne))
allocate(nq_neutral(Nq, Nq, Nq, Ne))

nenq = func2quad(nodes, elems, xiq, fnen)
nq_pos = func2quad(nodes, elems, xiq, fn_pos)
! Make the charge density net neutral (zero integral):
background = integral(nodes, elems, wtq3, nq_pos) / (Lx*Ly*Lz)
print *, "Total (positive) electronic charge: ", background * (Lx*Ly*Lz)
print *, "Subtracting constant background (Q/V): ", background
nq_neutral = nq_pos - background
!call assemble_3d(xin, nodes, elems, ib, xiq, wtq3, phihq, dphihq, &
!    4*pi*nq_neutral, Ap, Aj, Ax, rhs)


Nesub = count(mask_elems == 1)

Nbelem = (p+1)**3
Asize = Nesub * Nbelem * (Nbelem+1)/2

allocate(matAi(Asize))
allocate(matAj(Asize))
allocate(matAx(Asize))

elx = nodes(1, elems(7, 1)) - nodes(1, elems(1, 1)) ! Element sizes
ely = nodes(2, elems(7, 1)) - nodes(2, elems(1, 1))
elz = nodes(3, elems(7, 1)) - nodes(3, elems(1, 1))
print *, "precalculating"
allocate(global_to_local(Nb))
print *, "assembling 1"
call assemble_3d_coo_A_subdomain_spectral(Ne, p, ib, dphihq, elx, ely, elz, &
    wtq, mask_elems, matAi, matAj, matAx, idx, jac_det, global_to_local)
call assemble_3d_coo_rhs_spectral(Ne, p, 4*pi*nq_neutral, jac_det, wtq3, &
    ib, rhs)
print *, "preparing data..."
Nbsub = count(global_to_local /= 0)
allocate(local_to_global(Nbsub))
j = 1
do i = 1, size(global_to_local)
    if (global_to_local(i) > 0) then
        local_to_global(j) = i
        j = j + 1
    end if
end do
allocate(elems_sub(Nesub*Nbelem))
j = 1
do i = 1, Ne
    if (mask_elems(i) == 1) then
        elems_sub(j:j+Nbelem-1) = global_to_local(reshape(ib(:, :, :, i), [Nbelem]))
        j = j + Nbelem
    end if
end do
allocate(isegns(Nesub))
j = 1
do i = 1, Ne
    if (mask_elems(i) == 1) then
        isegns(j) = i
        j = j + 1
    end if
end do
allocate(nnets(Nesub))
nnets = Nbelem
allocate(nndfs(Nbsub))
nndfs = 1
jacx = elx/2
jacy = ely/2
jacz = elz/2
allocate(xyz_global(3, Nb))
print *, "preparing xyz data"
do e = 1, Ne
    xp = (xin + 1) * jacx + nodes(1, elems(1, e))
    yp = (xin + 1) * jacy + nodes(2, elems(1, e))
    zp = (xin + 1) * jacz + nodes(3, elems(1, e))
    do iaz = 1, p+1
    do iay = 1, p+1
    do iax = 1, p+1
        i = ib(iax, iay, iaz, e)
        xyz_global(:, i) = [xp(iax), yp(iay), zp(iaz)]
    end do
    end do
    end do
end do
print *, "preparing the rest"
allocate(xyzs(Nbsub, 3))
xyzs = transpose(xyz_global(:, local_to_global))
allocate(ifixs(Nbsub))
ifixs = 0
allocate(fixvs(Nbsub))
fixvs = 0
allocate(rhss(Nbsub))
rhss = rhs(local_to_global)
allocate(sols(Nbsub))
sols = 0
matAi(:idx) = global_to_local(matAi(:idx))
matAj(:idx) = global_to_local(matAj(:idx))

nlevels = 2  ! bddc levels
nsublev = [nproc, 1]

! tell me how much subdomains should I load
nsub_loc_1 = -1
verbose_level = 1

matrixtype = 1

is_assembled_int = 0
is_rhs_complete_int = 1

print *, "bddcml_init"
call bddcml_init(nlevels, nsublev, size(nsublev), nsub_loc_1, comm, verbose_level, 1, 0)

isub = myid + 1

print *, "bddcml_upload_subdomain_data"
call bddcml_upload_subdomain_data(Ne, Nb, Nb, 3, 3, &
               isub, Nesub, Nbsub, Nbsub, &
               elems_sub,size(elems_sub),nnets,size(nnets), nndfs,size(nndfs), &
               local_to_global,size(local_to_global), local_to_global, &
                    size(local_to_global), isegns,size(isegns), &
               xyzs,size(xyzs, 1), size(xyzs, 2), &
               ifixs,size(ifixs), fixvs,size(fixvs), &
               rhss,size(rhss), is_rhs_complete_int, &
               sols,size(sols), &
               matrixtype, matAi, matAj, matAx, idx, is_assembled_int, &
               user_constraints, 0, 0, &
               element_data, 0, 0, &
               dof_data, 0)

print *, "bddcml_setup_preconditioner"
call bddcml_setup_preconditioner(matrixtype,&
                                   use_preconditioner_defaults, &
                                   parallel_division,&
                                   use_arithmetic_constraints,&
                                   use_adaptive_constraints,&
                                   use_user_constraints,&
                                   weights_type)

print *, "sum(rhs):    ", sum(rhs)
print *, "integral rhs:", integral(nodes, elems, wtq3, nq_neutral)
print *, "Solving..."

call bddcml_solve(comm, krylov_method, tol,maxit,ndecrmax, &
    recycling_int, max_number_of_stored_vectors, &
    num_iter, converged_reason, condition_number)

print *, "getting the solution"

call bddcml_download_global_solution(sol, size(sol))

print *, "broadcasting"

call MPI_BCAST(sol, size(sol), MPI_DOUBLE_PRECISION, 0, comm, ierr)

!print *, "myid = ", myid, "mask_elems = ", mask_elems
print *, "myid = ", myid, "Nbsub = ", Nbsub
!print *, "myid = ", myid, "local_to_global = ", local_to_global
!call MPI_BARRIER(comm, ierr)
!stop "OK"
!sol = solve(A, rhs)
!sol = solve_cg(Ap, Aj, Ax, rhs, zeros(size(rhs)), 1e-12_dp, 400, .true.)
print *, "fullc"
call c2fullc_3d(in, ib, sol, fullsol)
print *, "converting to quadrature points"
call fe2quad_3d_lobatto(elems, xin, in, fullsol, Vhq)

print *, "background"
background = integral(nodes, elems, wtq3, nenq) / (Lx*Ly*Lz)
print *, "Total (negative) ionic charge: ", background * (Lx*Ly*Lz)
print *, "Subtracting constant background (Q/V): ", background
nenq = nenq - background
!call assemble_3d(xin, nodes, elems, ib, xiq, wtq3, phihq, dphihq, &
!    4*pi*nenq, Ap, Aj, Ax, rhs)

print *, "assembly rhs"
call assemble_3d_coo_rhs_spectral(Ne, p, 4*pi*nenq, jac_det, wtq3, ib, rhs)
rhss = rhs(local_to_global)

print *, "update rhs"
call bddcml_change_subdomain_data(isub, &
    ifixs, size(ifixs), fixvs, size(fixvs), &
    rhss, size(rhss), is_rhs_complete_int, &
    sols, size(sols))

print *, "setup_data"
call bddcml_setup_new_data()

print *, "sum(rhs):    ", sum(rhs)
print *, "integral rhs:", integral(nodes, elems, wtq3, nenq)
print *, "Solving..."

call bddcml_solve(comm, krylov_method, tol,maxit,ndecrmax, &
    recycling_int, max_number_of_stored_vectors, &
    num_iter, converged_reason, condition_number)

print *, "getting the solution"
call bddcml_download_global_solution(sol, size(sol))

print *, "broadcasting"
call MPI_BCAST(sol, size(sol), MPI_DOUBLE_PRECISION, 0, comm, ierr)

print *, "done"
call bddcml_finalize()


!sol = solve_cg(Ap, Aj, Ax, rhs, zeros(size(rhs)), 1e-12_dp, 400)

print *, "fullc"

call c2fullc_3d(in, ib, sol, fullsol)
print *, "converting to quad points"
call fe2quad_3d_lobatto(elems, xin, in, fullsol, Venq)
!print *, "Saving Ven to VTK"
!call vtk_save("Venq.vtk", Nex, Ney, Nez, nodes, elems, xiq, Venq)
!print *, "Saving values of Ven on a line"
!call line_save("Venq_line.txt", xin, nodes, elems, in, fullsol, &
!    [0._dp, 0._dp, 0._dp], [1._dp, 0._dp, 0._dp], 500)
!print *, "    Done."

print *, "integrals..."

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
print *, "done"
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

real(dp) elemental function f(y)
! Function f(y) from Appendix A in [1].
!
! [1] Perrot, F. (1979). Gradient correction to the statistical electronic free
! energy at nonzero temperatures: Application to equation-of-state
! calculations. Physical Review A, 20(2), 586–594.
real(dp), intent(in) :: y ! must be positive
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
if (y <= y0) then
    f = log(y)
    do i = 0, 7
        f = f + c(i+1) * y**i
    end do
else
    u = y**(2._dp / 3)
    f = d(1)*u
    do i = 1, 8
        f = f + d(i+1) / u**(2*i-1)
    end do
    ! Note: Few terms in [1] have "y" instead of "u" in them for y > y0, but
    ! that is obviously a typo.
end if
end function

end module


! ------------------------------------------------------------------------


program test_free_energy_bddc
use types, only: dp
use test_free_energy_bddc_utils, only: free_energy, read_pseudo
use constants, only: Ha2eV, pi
use utils, only: loadtxt, assert
use splines, only: spline3pars, iixmin, poly3
use interp3d, only: trilinear
use mpi_interface, only: mpi_finalize, mpi_comm_world, mpi_comm_rank, &
    mpi_comm_size, mpi_init
implicit none

real(dp) :: Eh, Een, Ts, Exc, Etot
integer :: p, DOF
real(dp) :: Z
real(dp) :: Rcut, L, T_eV, T_au
!  parallel variables
integer :: myid, comm_all, nproc, ierr


call MPI_INIT(ierr)
comm_all  = MPI_COMM_WORLD
call MPI_COMM_RANK(comm_all,myid,ierr)
call MPI_COMM_SIZE(comm_all,nproc,ierr)

!print *, "myid =", myid


Z = 1
Rcut = 0.3_dp
p = 6
L = 2
T_eV = 0.0862_dp
T_au = T_ev / Ha2eV
call free_energy(comm_all, L, 4, 4, 4, 2, 1, 1, p, T_au, nen, ne, Eh, Een, Ts, Exc, DOF)
Etot = Ts + Een + Eh + Exc
if (myid == 0) then
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

    ! The reference answers are converged to at least 1e-5:
    print *, abs(Ts  - (+10.61905))
    call assert(abs(Ts  - (+10.61905)) < 1e-4)
    print *, abs(Een - (- 3.80769))
    call assert(abs(Een - (- 3.80769)) < 1e-4)
    print *, abs(Eh  - (+ 1.30109))
    call assert(abs(Eh  - (+ 1.30109)) < 1e-4)
    print *, abs(Exc  - (- 1.43806))
    call assert(abs(Exc  - (- 1.43806)) < 1e-4)
end if

call MPI_FINALIZE(ierr)

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

real(dp) function ne(x, y, z) result(n)
real(dp), intent(in) :: x, y, z
real(dp), parameter :: alpha = 5, Z_ = 1
real(dp) :: r
r = sqrt(x**2+y**2+z**2)
n = Z_*alpha**3/pi**(3._dp/2)*exp(-alpha**2*R**2)
end function

end program
