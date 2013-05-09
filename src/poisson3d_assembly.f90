module poisson3d_assembly
use types, only: dp
use linalg, only: inv
use utils, only: assert, stop_error
use constants, only: pi
use sparse, only: coo2csr_canonical
implicit none
private
public assemble_3d, integral, func2quad, func_xyz

interface
    real(dp) function func_xyz(x, y, z)
    import :: dp
    implicit none
    real(dp), intent(in) :: x, y, z
    end function
end interface

contains

function func2quad(nodes, elems, xiq, func) result(fq)
! Return an array of function 'func' values at quadrature points
real(dp), intent(in):: nodes(:, :), xiq(:)
procedure(func_xyz) :: func
integer, intent(in):: elems(:, :)
integer :: Ne, e, iqx, iqy, iqz
real(dp), dimension(size(xiq), size(xiq), size(xiq), size(elems, 2)) :: fq
real(dp), dimension(size(xiq)) :: x, y, z, xp, yp, zp
real(dp) :: lx, ly, lz
real(dp) :: jacx, jacy, jacz, jac_det
lx = nodes(1, elems(7, 1)) - nodes(1, elems(1, 1)) ! Element sizes
ly = nodes(2, elems(7, 1)) - nodes(2, elems(1, 1))
lz = nodes(3, elems(7, 1)) - nodes(3, elems(1, 1))
jacx = lx/2
jacy = ly/2
jacz = lz/2
jac_det = abs(jacx*jacy*jacz)
xp = (xiq + 1) * jacx
yp = (xiq + 1) * jacy
zp = (xiq + 1) * jacz
Ne = size(elems, 2)
do e = 1, Ne
    x = xp + nodes(1, elems(1, e))
    y = yp + nodes(2, elems(1, e))
    z = zp + nodes(3, elems(1, e))
    do iqz = 1, size(xiq)
    do iqy = 1, size(xiq)
    do iqx = 1, size(xiq)
        fq(iqx, iqy, iqz, e) = func(x(iqx), y(iqy), z(iqz))
    end do
    end do
    end do
end do
end function

subroutine assemble_3d(xin, nodes, elems, ib, xiq, wtq, phihq, dphihq, &
        rhsq, matBp, matBj, matBx, rhs)
! Assemble Poisson equation on a 3D hexahedral uniform mesh
! It solves:
!   -\Nabla^2 V(x, y, z) = f(x, y, z)
! where f(x, y, z) is specified in "rhsq". Typically, one would use f = 4*pi*n,
! where "n" is the positive particle density.
real(dp), intent(in):: xin(:), nodes(:, :), xiq(:), wtq(:, :, :), &
    phihq(:, :), dphihq(:, :), rhsq(:, :, :, :)
integer, intent(in):: elems(:, :), ib(:, :, :, :)
real(dp), intent(out):: rhs(:)
integer, allocatable, intent(out) :: matBp(:), matBj(:)
real(dp), allocatable, intent(out) :: matBx(:)
integer, allocatable :: matAi(:), matAj(:)
real(dp), allocatable :: matAx(:)
integer :: Ne, p, e, i, j, iqx, iqy, iqz
real(dp), dimension(size(xiq), size(xiq), size(xiq), &
    size(xin), size(xin), size(xin)) :: phi_v, phi_dx, phi_dy, phi_dz
real(dp), dimension(size(xiq), size(xiq), size(xiq)) :: fq
real(dp) :: lx, ly, lz
integer :: ax, ay, az, bx, by, bz
real(dp) :: jacx, jacy, jacz, jac_det
real(dp), dimension(size(xiq), size(xiq), size(xiq), &
    size(xiq), size(xiq), size(xiq)) :: Am_loc
integer :: idx, maxidx

Ne = size(elems, 2)
p = size(xin) - 1
! 3D shape functions
do az = 1, p+1
do ay = 1, p+1
do ax = 1, p+1
    do iqz = 1, size(xiq)
    do iqy = 1, size(xiq)
    do iqx = 1, size(xiq)
        phi_v (iqx, iqy, iqz, ax, ay, az) = &
             phihq(iqx, ax) *  phihq(iqy, ay) *  phihq(iqz, az)
        phi_dx(iqx, iqy, iqz, ax, ay, az) = &
            dphihq(iqx, ax) *  phihq(iqy, ay) *  phihq(iqz, az)
        phi_dy(iqx, iqy, iqz, ax, ay, az) = &
             phihq(iqx, ax) * dphihq(iqy, ay) *  phihq(iqz, az)
        phi_dz(iqx, iqy, iqz, ax, ay, az) = &
             phihq(iqx, ax) *  phihq(iqy, ay) * dphihq(iqz, az)
    end do
    end do
    end do
end do
end do
end do
rhs=0
! Precalculate as much as possible:
lx = nodes(1, elems(7, 1)) - nodes(1, elems(1, 1)) ! Element sizes
ly = nodes(2, elems(7, 1)) - nodes(2, elems(1, 1))
lz = nodes(3, elems(7, 1)) - nodes(3, elems(1, 1))
jacx = lx/2
jacy = ly/2
jacz = lz/2
jac_det = abs(jacx*jacy*jacz)
phi_dx = phi_dx / jacx
phi_dy = phi_dy / jacy
phi_dz = phi_dz / jacz
! Precalculate element matrices:
print *, "Precalculating local element matrix..."
do bz = 1, p+1
do by = 1, p+1
do bx = 1, p+1
    do az = 1, p+1
    do ay = 1, p+1
    do ax = 1, p+1
        Am_loc(ax, ay, az, bx, by, bz) = sum(( &
            phi_dx(:, :, :, ax, ay, az)*phi_dx(:, :, :, bx, by, bz) + &
            phi_dy(:, :, :, ax, ay, az)*phi_dy(:, :, :, bx, by, bz) + &
            phi_dz(:, :, :, ax, ay, az)*phi_dz(:, :, :, bx, by, bz)) &
            * jac_det * wtq)
    end do
    end do
    end do
end do
end do
end do
print *, "Assembly..."
idx = 0
maxidx = Ne*(p+1)**6
print *, "Number of COO matrix entries:", maxidx
allocate(matAi(maxidx), matAj(maxidx), matAx(maxidx))
do e = 1, Ne
    fq = rhsq(:, :, :, e)
    fq = fq * jac_det * wtq
    do bz = 1, p+1
    do by = 1, p+1
    do bx = 1, p+1
        j = ib(bx, by, bz, e)
        if (j==0) cycle
        do az = 1, p+1
        do ay = 1, p+1
        do ax = 1, p+1
            i = ib(ax, ay, az, e)
            if (i == 0) cycle
            if (j > i) cycle
            idx = idx + 1
            matAi(idx) = i
            matAj(idx) = j
            matAx(idx) = Am_loc(ax, ay, az, bx, by, bz)
            if (i /= j) then
                ! Symmetric contribution
                idx = idx + 1
                matAi(idx) = j
                matAj(idx) = i
                matAx(idx) = Am_loc(ax, ay, az, bx, by, bz)
            end if
        end do
        end do
        end do
        rhs(j) = rhs(j) + sum(phi_v(:, :, :, bx, by, bz) * fq)
    end do
    end do
    end do
end do
print *, "Converting COO -> CSR..."
call coo2csr_canonical(matAi(:idx), matAj(:idx), matAx(:idx), &
    matBp, matBj, matBx)
print *, "CSR matrix, dimension:", size(matBp)-1
print *, "CSR matrix, number of nonzeros:", size(matBx), &
    "density", real(size(matBx), dp) / (size(matBp)-1)**2
end subroutine

real(dp) function integral(nodes, elems, wtq, fq) result(r)
real(dp), intent(in) :: nodes(:, :)
integer, intent(in) :: elems(:, :)
real(dp), intent(in) :: wtq(:, :, :), fq(:, :, :, :)
real(dp) :: jacx, jacy, jacz, l(3), jac_det
integer :: e, Ne
Ne = size(elems, 2)
r = 0
do e = 1, Ne
    ! l is the diagonal vector:
    l = nodes(:, elems(7, e)) - nodes(:, elems(1, e))
    ! Assume rectangular shape:
    jacx = l(1)/2
    jacy = l(2)/2
    jacz = l(3)/2
    jac_det = abs(jacx*jacy*jacz)
    r = r + sum(fq(:, :, :, e) * jac_det * wtq)
end do
end function

end module

