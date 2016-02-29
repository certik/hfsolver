module ewald_sums
use types, only: dp
use constants, only: pi
implicit none
private
public ewald, ewald2, direct_sum, fred2fcart, sred2scart, ewald_box, &
    fft_neutralized

contains

subroutine ewald(eew,gmet,grewtn,natom,ntypat,rmet,typat,ucvol,xred,zion)
! This subroutine was taken from ABINIT
!
!{\src2tex{textfont=tt}}
!!****f* ABINIT/ewald
!!
!! NAME
!! ewald
!!
!! FUNCTION
!! Compute Ewald energy and derivatives with respect to dimensionless
!!  reduced atom coordinates xred.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gmet(3,3)=metric tensor in reciprocal space (bohr^-2)
!! natom=number of atoms in unit cell
!! ntypat=numbe of type of atoms
!! rmet(3,3)=metric tensor in real space (bohr^2)
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume (bohr^3)
!! xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!! zion(ntypat)=charge on each type of atom (real number)
!!
!! OUTPUT
!! eew=final ewald energy in hartrees
!! grewtn(3,natom)=grads of eew wrt xred(3,natom), hartrees.
!Arguments ------------------------------------
!scalars
integer, intent(in) :: natom, ntypat
real(dp), intent(in) :: ucvol
real(dp), intent(out) :: eew
!arrays
integer, intent(in) :: typat(natom)
real(dp), intent(in) :: gmet(3,3), rmet(3,3), xred(3,natom), zion(ntypat)
real(dp), intent(out) :: grewtn(3, natom)

!Local variables-------------------------------
integer :: ia,ib,ig1,ig2,ig3,ir1,ir2,ir3,newg,newr,ng,nr
real(dp) :: arg,c1i,ch,chsq,drdta(3),eta,fac
real(dp) :: fraca(3),fracb(3),gsq,gsum,phi,phr,r(3)
real(dp) :: reta,rmagn,rsq,sumg,summi,summr,sumr
real(dp) :: term

!Add up total charge and sum of $charge^2$ in cell
ch   = sum(zion(typat))
chsq = sum(zion(typat)**2)

!Compute eta, the Ewald summation convergence parameter,
!for approximately optimized summations:
!A bias is introduced, because G-space summation scales
!better than r space summation ! Note : debugging is the most
!easier at fixed eta.
eta = pi*200._dp/33._dp*sqrt(1.69_dp*sum(gmet)/sum(rmet))

!Conduct reciprocal space summations
fac = pi**2/eta
gsum = 0
grewtn = 0
!Sum over G space, done shell after shell until all
!contributions are too small.
ng = 0
do
    ng = ng + 1
    newg = 0
    do ig3 = -ng, ng
    do ig2 = -ng, ng
    do ig1 = -ng, ng
        ! Exclude shells previously summed over
        if (abs(ig1)==ng .or. abs(ig2)==ng .or. abs(ig3)==ng .or. ng==1) then
            ! gsq is G dot G = |G|^2
            gsq = gmet(1,1)*ig1**2+gmet(2,2)*ig2**2+gmet(3,3)*ig3**2 + &
                2*(gmet(2,1)*ig1*ig2+gmet(3,1)*ig1*ig3+gmet(3,2)*ig3*ig2)
            ! Skip g=0:
            if (gsq > 1e-20_dp) then
                arg = fac*gsq
                ! Larger arg gives 0 contribution because of exp(-arg)
                if (arg <= 80) then
                    ! When any term contributes then include next shell
                    newg = 1
                    term = exp(-arg)/gsq
                    summr = 0
                    summi = 0
                    ! Note that if reduced atomic coordinates xred drift
                    ! outside of unit cell (outside [0,1)) it is irrelevant in
                    ! the following term, which only computes a phase.
                    do ia = 1, natom
                        arg = 2*pi*dot_product([ig1, ig2, ig3], xred(:, ia))
                        ! Sum real and imaginary parts (avoid complex variables)
                        summr = summr+zion(typat(ia))*cos(arg)
                        summi = summi+zion(typat(ia))*sin(arg)
                    end do
                    ! The following two checks avoid an annoying
                    ! underflow error message
                    if (abs(summr) < 1e-16_dp) summr = 0
                    if (abs(summi) < 1e-16_dp) summi = 0
                    ! The product of term and summr**2 or summi**2 below
                    ! can underflow if not for checks above
                    gsum = gsum + term*(summr**2+summi**2)
                    do ia = 1, natom
                        ! Again only phase is computed so xred may fall outside
                        ! [0,1).
                        arg = 2*pi*dot_product([ig1, ig2, ig3], xred(:, ia))
                        phr =  cos(arg)
                        phi = -sin(arg)
                        ! (note: do not need real part, commented out)
                        ! c1r=(phr*summr-phi*summi)*(term*zion(typat(ia)))
                        c1i = (phi*summr+phr*summi)*(term*zion(typat(ia)))
                        ! compute coordinate gradients
                        grewtn(:, ia) = grewtn(:, ia) - c1i*[ig1, ig2, ig3]
                    end do
                end if
            end if
        end if
    end do
    end do
    end do
    ! Check if new shell must be calculated
    if (newg == 0) exit
end do
sumg = gsum/(2*pi*ucvol)
!Stress tensor is now computed elsewhere (ewald2) hence do not need
!length scale gradients (used to compute them here).

!normalize coordinate gradients by unit cell volume ucvol
grewtn = -2*grewtn(:,:)/ucvol

!Conduct real space summations
reta = sqrt(eta)
fac = 2*sqrt(eta/pi)
sumr = 0
!In the following a summation is being conducted over all
!unit cells (ir1, ir2, ir3) so it is appropriate to map all
!reduced coordinates xred back into [0,1).
!
!Loop on shells in r-space as was done in g-space
nr = 0
do
    nr = nr + 1
    newr = 0
    ! Triple loop over real space points and associated condition of new shell
    do ir3 = -nr, nr
    do ir2 = -nr, nr
    do ir1 = -nr, nr
      if (abs(ir3)==nr .or. abs(ir2)==nr .or. abs(ir1)==nr .or. nr==1) then
        do ia = 1, natom
          ! Map reduced coordinate xred(mu,ia) into [0,1)
          fraca = xred(:,ia)-aint(xred(:,ia))+0.5_dp-sign(0.5_dp,xred(:,ia))
          drdta = 0
          do ib = 1, natom
              fracb = xred(:,ib)-aint(xred(:,ib))+0.5_dp-sign(0.5_dp,xred(:,ib))
              r = [ir1, ir2, ir3] + fracb - fraca
              rsq = rmet(1,1)*r(1)**2+rmet(2,2)*r(2)**2+rmet(3,3)*r(3)**2 + &
                2*(rmet(2,1)*r(2)*r(1)+rmet(3,2)*r(3)*r(2)+rmet(3,1)*r(1)*r(3))
              ! Avoid zero denominators in 'term':
              if (rsq >= 1e-24_dp) then
                  ! Note: erfc(8) is about 1.1e-29,
                  ! so do not bother with larger arg.
                  ! Also: exp(-64) is about 1.6e-28,
                  ! so do not bother with larger arg**2 in exp.
                  term=0
                  if (eta*rsq < 64) then
                      newr = 1
                      rmagn = sqrt(rsq)
                      term = erfc(reta*rmagn) / rmagn
                      sumr = sumr+zion(typat(ia))*zion(typat(ib))*term
                      term = zion(typat(ia)) * zion(typat(ib)) * &
                          (term+fac*exp(-eta*rsq))/rsq
                      ! Compute terms related to coordinate gradients
                      drdta = drdta + term*matmul(rmet, r)
                  end if
              end if
          end do ! ib
          grewtn(:, ia) = grewtn(:, ia) + drdta
        end do ! ia
      end if
    end do
    end do
    end do
    ! Check if new shell must be calculated
    if (newr == 0) exit
end do
! Finally assemble Ewald energy, eew
eew = sumg+sumr/2-chsq*reta/sqrt(pi)-pi*ch**2/(2*eta*ucvol)
end subroutine ewald


subroutine ewald2(gmet,natom,ntypat,rmet,rprimd,gprimd,stress,typat,ucvol, &
        xred,zion)
! This subroutine was taken from ABINIT
!
!{\src2tex{textfont=tt}}
!!****f* ABINIT/ewald2
!!
!! NAME
!! ewald2
!!
!! FUNCTION
!! Compute the part of the stress tensor coming from the Ewald energy
!! which is calculated by derivating the Ewald energy with respect to
!! strain.
!! See Nielsen and Martin, Phys. Rev. B 32, 3792 (1985).
!! Definition of stress tensor is $(1/ucvol)*d(Etot)/d(strain(a,b))$.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2014 ABINIT group (JCC, DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gmet(3,3)=metric tensor in reciprocal space (bohr^-2)
!! natom=number of atoms in umit cell
!! ntypat=number of type of atoms
!! rmet(3,3)=metric tensor in real space (bohr^2) (inverse transpose of gmet)
!! rprimd(3,3)=dimensional primitive translations in real space (bohr)
!! gprimd(3,3)=dimensional primitive translations in reciprocal space (1/bohr)
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume (bohr^3)
!! xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!! zion(ntypat)=charge on each type of atom (real number)
!!
!! OUTPUT
!! $stress(6)=(1/ucvol)*gradient$ of Ewald energy with respect to strain,
!!      in hartrees/bohr^3
!! Cartesian components of stress are provided for this symmetric
!! tensor in the order 11 22 33 32 31 21.
!!
!Arguments ------------------------------------
!scalars
integer, intent(in) :: natom, ntypat
real(dp), intent(in) :: ucvol
!arrays
integer, intent(in) :: typat(natom)
real(dp), intent(in) :: gmet(3,3), rmet(3,3), rprimd(3,3), gprimd(3,3), &
    xred(3,natom), zion(ntypat)
real(dp), intent(out) :: stress(6)

!Local variables-------------------------------
integer :: ia,ib,ig1,ig2,ig3,ir1,ir2,ir3,newg,newr,ng,nr
real(dp) :: arg1,arg2,arg3,ch,dderfc,eta,fac,fraca(3)
real(dp) :: fracb(3),g(3),gsq,r(3),rc(3)
real(dp) :: reta,rmagn,rsq,summi,summr,t(6),term1
real(dp) :: term4
real(dp) :: strg(6),strr(6)

!Add up total charge and sum of charge^2 in cell
ch = sum(zion(typat))

!Compute eta, the Ewald summation convergence parameter,
!for approximately optimized summations:
!Here, a bias is introduced, because G-space summation scales
!better than r space summation !
eta = pi*200.0_dp/33.0_dp*sqrt(1.69_dp*sum(gmet)/sum(rmet))

fac = pi**2/eta

!Conduct reciprocal space summations
strg = 0
!Sum over G space, done shell after shell until all
!contributions are too small
ng=0
do
    ng = ng+1
    newg = 0
    do ig3 = -ng, ng
    do ig2 = -ng, ng
    do ig1 = -ng, ng
        ! Exclude shells previously summed over
        if (abs(ig1)==ng .or. abs(ig2)==ng .or. abs(ig3)==ng .or. ng==1) then
            g = matmul(gprimd, [ig1, ig2, ig3]) ! Cartesian components of each G
            gsq = sum(g**2) ! |G|^2 (no pi factors)
            ! skip g=0:
            if (gsq > 1e-20_dp) then
                arg1 = fac*gsq
                ! larger arg1 gives 0 contribution because of exp(-arg1)
                if (arg1 <= 80) then
                    ! When any term contributes then include next shell
                    newg = 1
                    term1 = exp(-arg1)/arg1
                    summr = 0
                    summi = 0
                    do ia = 1, natom
                        arg2 = 2*pi*dot_product([ig1, ig2, ig3], xred(:, ia))
                        ! Sum real and imaginary parts (avoid complex variables)
                        summr = summr+zion(typat(ia))*cos(arg2)
                        summi = summi+zion(typat(ia))*sin(arg2)
                    end do
                    ! Avoid underflow error messages
                    if (abs(summr) < 1e-16_dp) summr = 0
                    if (abs(summi) < 1e-16_dp) summi = 0
                    t = (2/gsq)*(1+arg1) * [g(1)**2, g(2)**2, g(3)**2, &
                        g(2)*g(3), g(1)*g(3), g(1)*g(2)] - [1, 1, 1, 0, 0, 0]
                    strg = strg + t*term1*(summr*summr+summi*summi)
                end if
           end if
        end if
    end do
    end do
    end do
    ! Check if new shell must be calculated
    if (newg == 0) exit
end do

!Conduct real space summations
reta = sqrt(eta)
strr = 0
!Loop on shells in r-space as was done in g-space
nr = 0
do
    nr = nr + 1
    newr = 0
    ! Triple loop over real space points, and associated new shell condition
    do ir3 = -nr, nr
    do ir2 = -nr, nr
    do ir1 = -nr, nr
      if (abs(ir3)==nr .or. abs(ir2)==nr .or. abs(ir1)==nr .or. nr==1) then
        do ia=1,natom
          ! Convert reduced atomic coordinates to [0,1)
          fraca = xred(:,ia)-aint(xred(:,ia))+0.5_dp-sign(0.5_dp,xred(:,ia))
          do ib = 1, natom
              fracb=xred(:,ib)-aint(xred(:,ib))+0.5_dp-sign(0.5_dp,xred(:,ib))
              r = [ir1, ir2, ir3] + fracb - fraca
              ! Convert from reduced to cartesian coordinates
              rc = matmul(rprimd, r)
              ! Compute |r|^2
              rsq = sum(rc**2)
              rmagn = sqrt(rsq)
              ! Avoid zero denominators in 'term':
              if (rmagn >= 1e-12_dp) then
                  ! Note: erfc(8) is about 1.1e-29,
                  ! so do not bother with larger arg.
                  ! Also: exp(-64) is about 1.6e-28,
                  ! so do not bother with larger arg**2 in exp.
                  arg3=reta*rmagn
                  if (arg3 < 8) then
                      newr = 1
                      ! derivative of the erfc:
                      dderfc = (-2/sqrt(pi))*exp(-eta*rsq)
                      term4 = zion(typat(ia)) * zion(typat(ib)) * &
                          (dderfc-erfc(arg3)/arg3)
                      strr = strr + term4 * [rc(1)**2, rc(2)**2, rc(3)**2, &
                          rc(2)*rc(3), rc(1)*rc(3), rc(1)*rc(2)] / rsq
                  end if
              end if
          end do ! ib
        end do ! ia
      end if
    end do ! ir1
    end do ! ir2
    end do ! ir3
    ! Check if new shell must be calculated
    if (newr == 0) exit
end do
! Finally assemble stress tensor coming from Ewald energy, stress
! (note division by unit cell volume in accordance with definition
! found in Nielsen and Martin, Phys. Rev. B 32, 3792 (1985).)
strg(:3) = strg(:3) + ch**2
stress = (reta*strr/2 + pi*strg/(2*ucvol*eta)) / ucvol
end subroutine ewald2


subroutine fred2fcart(fcart, fred, gprim)
! Convert forces from rprim reduced coordinates to cartesian coordinates
real(dp), intent(in) :: gprim(3, 3) ! gprim(:, i) = G_i (reciprocal vectors)
! fred(:, i) = F_i (force on i-th atom in gprim coordinates)
real(dp), intent(in) :: fred(:, :)
! fcart(:, i) = F_i (force on i-th atom in cartesian coordinates)
real(dp), intent(out) :: fcart(:, :)
integer :: i, n
n = size(fred, 2)
do i = 1, n
    fcart(:, i) = gprim(:, 1)*fred(1, i) + gprim(:, 2)*fred(2, i) &
        + gprim(:, 3)*fred(3, i)
end do
end subroutine

subroutine sred2scart(scart, sred, gprim)
! Convert stress tensor from rprim reduced coordinates to cartesian coordinates
! Symmetric storage mode for a 3x3 tensor is a 6 element array with elements
! 11, 22, 33, 32, 31 and 21.
real(dp), intent(in) :: gprim(3, 3) ! gprim(:, i) = G_i (reciprocal vectors)
! sred(6) stress tensor in reduced coordinates
real(dp), intent(in) :: sred(:)
! scart(6) stress tensor in cartesian coordinates
real(dp), intent(out) :: scart(:)
real(dp) :: w(3, 3)
w(1, 1) = sred(1)
w(2, 2) = sred(2)
w(3, 3) = sred(3)
w(3, 2) = sred(4); w(2, 3) = sred(4)
w(3, 1) = sred(5); w(1, 3) = sred(5)
w(2, 1) = sred(6); w(1, 2) = sred(6)
w = matmul(matmul(gprim, w), transpose(gprim))
scart(1) = w(1, 1)
scart(2) = w(2, 2)
scart(3) = w(3, 3)
scart(4) = w(2, 3)
scart(5) = w(1, 3)
scart(6) = w(1, 2)
end subroutine

subroutine direct_sum(q, r, L, ncut, E, forces)
! Calculates the Coulombic energy E and forces as a direc lattice sum.
!
! We use a spherically ordered direct lattice sum. If the unit cell dipole
! moment is non-zero, the sum needs to be corrected in order to obtain the
! correct Ewald energy and forces, see [1] for more details.
!
! The summation is converging very slowly, but the subroutine is simple to
! implement and maintain, so it is used for testing the correctness of faster
! more advanced methods. The convergence is achieved by increasing 'ncut'.
!
! [1] Roberts, J. E., Schnitker, J. (1994). How the unit cell surface charge
! distribution affects the energetics of ion–solvent interactions in
! simulations. The Journal of Chemical Physics, 101(6), 5024.
! doi:10.1063/1.467425
real(dp), intent(in) :: q(:) ! q(i) is charge of i-th particle
real(dp), intent(in) :: r(:, :) ! r(:, i) is (x, y, z) coordinates of i-th par.
real(dp), intent(in) :: L ! length of the unit cell (box)
integer, intent(in) :: ncut ! cutoff
real(dp), intent(out) :: E ! Calculated energy
real(dp), intent(out) :: forces(:, :) ! forces(:, i) is a force on i-th particle
integer :: N, i, j, nx, ny, nz
real(dp) :: d_ji(3), d
N = size(q)
E = 0
forces = 0

! Calculate the spherically ordered direct lattice sum
do i = 1, N
    do j = 1, N
        do nx = -ncut, ncut
        do ny = -ncut, ncut
        do nz = -ncut, ncut
            if (nx == 0 .and. ny == 0 .and. nz == 0 .and. i == j) cycle
            if (sqrt(real(nx**2+ny**2+nz**2, dp)) > ncut) cycle
            ! Vector pointing from particle j -> i: d_ji = r_i - r_j
            d_ji = r(:, i)-r(:, j) + [nx, ny, nz]*L
            d = sqrt(sum(d_ji**2)) ! |r_i - r_j|
            E = E + q(i)*q(j) / d
            forces(:, i) = forces(:, i) + q(i)*q(j)*d_ji/d**3
        end do
        end do
        end do
    end do
end do
E = E / 2

! Calculate the dipole moment of the unit cell
d_ji = 0
do i = 1, N
    d_ji = d_ji + q(i)*r(:, i)
end do

! Correct the sum in order to obtain Ewald energy and forces
! Eq. (9) in [1]
E = E - 2*pi/3*sum(d_ji**2)/L**3
! Eq. (13) in [1]
do i = 1, N
    forces(:, i) = forces(:, i) + 4*pi/(3*L**3)*q(i)*d_ji
end do
end subroutine

subroutine ewald_box(L, x, q, E, forces, stress)
! Special case of Ewald summation for a box L^3. This calls ewald() and
! ewald2() internally. Because of the simplified geometry, it has simpler
! arguments and it calculates the various geometric quantities automatically
! before calling ewald(). All quantities are in a.u.
real(dp), intent(in) :: L ! Length of the box
real(dp), intent(in) :: x(:, :) ! x(:, i) position of i-th ion in [0, L]^3
real(dp), intent(in) :: q(:) ! r(i) charge of i-th ion
real(dp), intent(out) :: E ! ion-ion electrostatic potential energy
real(dp), intent(out) :: forces(:, :) ! forces(:, i) forces on i-th ion
real(dp), intent(out) :: stress(:) ! stress(6) the stress tensor

real(dp) :: rmet(3, 3), gmet(3, 3), rprim(3, 3), gprim(3, 3), ucvol
real(dp), dimension(3, size(x, 2)) :: xred, grewtn
integer :: typat(size(x, 2))
real(dp) :: zion(size(x, 2))
integer :: natom, ntypat, i
natom = size(x, 2)
ntypat = natom
do i = 1, ntypat
    typat(i) = i
end do
zion = q

! Setup geometric quantities:
rmet = 0
rmet(1, 1) = L**2
rmet(2, 2) = L**2
rmet(3, 3) = L**2

! gmet = inv(rmet)
gmet = 0
gmet(1, 1) = 1/L**2
gmet(2, 2) = 1/L**2
gmet(3, 3) = 1/L**2

! ucvol = sqrt(det(rmet))
ucvol = L**3

! Reciprocal primitive vectors (without 2*pi) in cartesian coordinates.
! gmet = matmul(transpose(gprim), gprim) = inv(rprim)
gprim = 0
gprim(1, 1) = 1 / L
gprim(2, 2) = 1 / L
gprim(3, 3) = 1 / L

! Real space primitive vectors
! rmet = matmul(transpose(rprim), rprim)
rprim = 0
rprim(1, 1) = L
rprim(2, 2) = L
rprim(3, 3) = L

xred = x / L

call ewald(E,gmet,grewtn,natom,ntypat,rmet,typat,ucvol,xred,zion)
call ewald2(gmet,natom,ntypat,rmet,rprim,gprim,stress,typat,ucvol,xred,zion)
call fred2fcart(forces, grewtn, gprim)
! Make forces out of gradients:
forces = -forces
stress = -stress * ucvol
end subroutine

subroutine fft_neutralized(L, x, q, E)
use ofdft_fft, only: reciprocal_space_vectors, real2fourier, &
    real_space_vectors, fourier2real
real(dp), intent(in) :: L ! Length of the box
real(dp), intent(in) :: x(:, :) ! x(:, i) position of i-th ion in [0, L]^3
real(dp), intent(in) :: q(:) ! r(i) charge of i-th ion
real(dp), intent(out) :: E ! ion-ion electrostatic potential energy
integer :: N, ii, i, j, k, Ng
real(dp) :: rc, Ig, v0, Isph, rho_minus, r, Xj(3), d(3)
real(dp), allocatable :: rho_tilde_minus(:, :, :)
complex(dp), allocatable :: rho_tilde_minusG(:, :, :)
real(dp), allocatable :: G(:, :, :, :), G2(:, :, :), Xn(:,:,:,:)

rho_minus = -sum(q)/L**3
!rc = L/2
rc = 0.11_dp
Ng = 256

!Ig = 10976 / (17875*rc)
!Isph = 14*pi*rc**2/75
!v0 = 12/(5*rc)

Ig = (517399110137._dp/809268832200._dp)/rc
Isph = 91*pi*rc**2/690
v0 = 11/(4*rc)


N = size(q)
E = 0

do i = 1, N
    E = E + 1/2._dp * q(i)**2 * (Ig-v0) + q(i) * Isph * rho_minus
end do

print *, E

allocate(G(Ng, Ng, Ng, 3), G2(Ng, Ng, Ng), Xn(Ng, Ng, Ng, 3))
allocate(rho_tilde_minus(Ng,Ng,Ng), rho_tilde_minusG(Ng,Ng,Ng))

call real_space_vectors(L, Xn)
call reciprocal_space_vectors(L, G, G2)

rho_tilde_minus = rho_minus
do ii = 1, N
    do k = 1, Ng
    do j = 1, Ng
    do i = 1, Ng
        Xj = x(:, ii)-Xn(i,j,k,:)+[L/2, L/2, L/2]
        Xj = Xj - L*floor(Xj/L)
        d = [L/2, L/2, L/2] - Xj
        r = sqrt(sum(d**2))
        rho_tilde_minus(i,j,k) = rho_tilde_minus(i,j,k) &
            + q(ii)*g3_fn(r, rc)
    end do
    end do
    end do
end do

call real2fourier(rho_tilde_minus, rho_tilde_minusG)
!rho_tilde_minusG = rho_tilde_minusG*4*pi/G2
!rho_tilde_minusG(1,1,1) = 0
!call fourier2real(rho_tilde_minusG, rho_tilde_minus)

!print *, "charge density along diagonal"
!print *, "potential along diagonal"
!do i = 1, Ng
!    print *, i, sqrt(sum(Xn(i,i,i,:)**2)), rho_tilde_minus(i,i,i)
!end do

do k = 1, Ng
do j = 1, Ng
do i = 1, Ng
    if (i == 1 .and. j == 1 .and. k == 1) cycle
    E = E + 2*pi*abs(rho_tilde_minusG(i,j,k))**2/G2(i,j,k)*L**3
end do
end do
end do

print *, E

contains

    real(dp) function g_fn(r, rc) result(g)
    real(dp), intent(in) :: r, rc
    if (r > rc) then
        g = 0
    else
        g = -21*(r-rc)**3*(6*r**2+3*r*rc+rc**2)/(5*pi*rc**8)
    end if
    end function

    real(dp) function g2_fn(r, rc) result(g)
    real(dp), intent(in) :: r, rc
    real(dp), parameter :: C = 0.4410888872766044004562838172_dp
    if (r >= rc) then
        g = 0
    else
        ! Bump function
        g = exp(-1/(1-(r/rc)**2)) / (C*rc**3)
    end if
    end function

    real(dp) function g3_fn(r, rc) result(g)
    real(dp), intent(in) :: r, rc
    if (r >= rc) then
        g = 0
    else
        g = 21*(r - rc)**10*(48620*r**9 + 24310*r**8*rc + 11440*r**7*rc**2 + &
        5005*r**6*rc**3 + 2002*r**5*rc**4 + 715*r**4*rc**5 + 220*r**3*rc**6 + &
        55*r**2*rc**7 + 10*r*rc**8 + rc**9)/(4*pi*rc**22)
    end if
    end function

end subroutine

end module
