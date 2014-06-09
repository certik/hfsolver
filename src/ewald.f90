module ewald_sums
use types, only: dp
use constants, only: pi
implicit none
private
public ewald, ewald2, direct_sum, fred2fcart, sred2scart, ewald_box

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
!!
!! PARENTS
!!      setvtr
!!
!! CHILDREN
!!      derfc,wrtout
!!
!! SOURCE


!Arguments ------------------------------------
!scalars
 real(dp), parameter :: two_pi = 2*pi
 integer,intent(in) :: natom,ntypat
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: eew
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: gmet(3,3),rmet(3,3),xred(3,natom),zion(ntypat)
 real(dp),intent(out) :: grewtn(3,natom)

!Local variables-------------------------------
!scalars
 integer :: ia,ib,ig1,ig2,ig3,ir1,ir2,ir3,newg,newr,ng,nr
 real(dp) :: arg,c1i,ch,chsq,derfc_arg,direct,drdta1,drdta2,drdta3,eta,fac
 real(dp) :: fraca1,fraca2,fraca3,fracb1,fracb2,fracb3,gsq,gsum,phi,phr,r1
 real(dp) :: r1a1d,r2,r2a2d,r3,r3a3d,recip,reta,rmagn,rsq,sumg,summi,summr,sumr
 real(dp) :: t1,term
! character(len=500) :: message

! *************************************************************************

!DEBUG
!write(std_out,*)' ewald : enter '
!stop
!ENDDEBUG

!Add up total charge and sum of $charge^2$ in cell
 chsq=0._dp
 ch=0._dp
 do ia=1,natom
   ch=ch+zion(typat(ia))
   chsq=chsq+zion(typat(ia))**2
 end do

!Compute eta, the Ewald summation convergence parameter,
!for approximately optimized summations:
 direct=rmet(1,1)+rmet(1,2)+rmet(1,3)+rmet(2,1)+&
& rmet(2,2)+rmet(2,3)+rmet(3,1)+rmet(3,2)+rmet(3,3)
 recip=gmet(1,1)+gmet(1,2)+gmet(1,3)+gmet(2,1)+&
& gmet(2,2)+gmet(2,3)+gmet(3,1)+gmet(3,2)+gmet(3,3)
!A bias is introduced, because G-space summation scales
!better than r space summation ! Note : debugging is the most
!easier at fixed eta.
 eta=pi*200.0_dp/33.0_dp*sqrt(1.69_dp*recip/direct)

!Conduct reciprocal space summations
 fac=pi**2/eta
 gsum=0._dp
 grewtn(:,:)=0.0_dp

!Sum over G space, done shell after shell until all
!contributions are too small.
 ng=0
 do
   ng=ng+1
   newg=0

   do ig3=-ng,ng
     do ig2=-ng,ng
       do ig1=-ng,ng

!        Exclude shells previously summed over
         if(abs(ig1)==ng .or. abs(ig2)==ng .or. abs(ig3)==ng&
&         .or. ng==1 ) then

!          gsq is G dot G = |G|^2
           gsq=gmet(1,1)*dble(ig1*ig1)+gmet(2,2)*dble(ig2*ig2)+&
&           gmet(3,3)*dble(ig3*ig3)+2._dp*(gmet(2,1)*dble(ig1*ig2)+&
&           gmet(3,1)*dble(ig1*ig3)+gmet(3,2)*dble(ig3*ig2))

!          Skip g=0:
           if (gsq>1.0d-20) then
             arg=fac*gsq

!            Larger arg gives 0 contribution because of exp(-arg)
             if (arg <= 80._dp) then
!              When any term contributes then include next shell
               newg=1
               term=exp(-arg)/gsq
               summr = 0.0_dp
               summi = 0.0_dp
!              Note that if reduced atomic coordinates xred drift outside
!              of unit cell (outside [0,1)) it is irrelevant in the following
!              term, which only computes a phase.
!              OCL SCALAR ! by MM for Fujitsu
               do ia=1,natom
                 arg=two_pi*(ig1*xred(1,ia)+ig2*xred(2,ia)+ig3*xred(3,ia))
!                Sum real and imaginary parts (avoid complex variables)
                 summr=summr+zion(typat(ia))*cos(arg)
                 summi=summi+zion(typat(ia))*sin(arg)
               end do

!              The following two checks avoid an annoying
!              underflow error message
               if (abs(summr)<1.d-16) summr=0.0_dp
               if (abs(summi)<1.d-16) summi=0.0_dp

!              The product of term and summr**2 or summi**2 below
!              can underflow if not for checks above
               t1=term*(summr*summr+summi*summi)
               gsum=gsum+t1

!              OCL SCALAR ! by MM for Fujitsu
               do ia=1,natom
!                Again only phase is computed so xred may fall outside [0,1).
                 arg=two_pi*(ig1*xred(1,ia)+ig2*xred(2,ia)+ig3*xred(3,ia))
                 phr= cos(arg)
                 phi=-sin(arg)
!                (note: do not need real part, commented out)
!                c1r=(phr*summr-phi*summi)*(term*zion(typat(ia)))
                 c1i=(phi*summr+phr*summi)*(term*zion(typat(ia)))
!                compute coordinate gradients
                 grewtn(1,ia)=grewtn(1,ia)-c1i*ig1
                 grewtn(2,ia)=grewtn(2,ia)-c1i*ig2
                 grewtn(3,ia)=grewtn(3,ia)-c1i*ig3
               end do

!              End condition of not larger than 80.0
             end if

!            End skip g=0
           end if

!          End triple loop over G s and associated new shell condition
         end if
       end do
     end do
   end do

!  Check if new shell must be calculated
   if (newg==0) exit

!  End the loop on ng (new shells). Note that there is one exit
!  from this loop.
 end do
!
 sumg=gsum/(two_pi*ucvol)

!Stress tensor is now computed elsewhere (ewald2) hence do not need
!length scale gradients (used to compute them here).

!normalize coordinate gradients by unit cell volume ucvol
 term=-2._dp/ucvol
 grewtn(:,:)=grewtn(:,:)*term
!call DSCAL(3*natom,term,grewtn,1)

!Conduct real space summations
 reta=sqrt(eta)
 fac=2._dp*sqrt(eta/pi)
 sumr=0.0_dp

!In the following a summation is being conducted over all
!unit cells (ir1, ir2, ir3) so it is appropriate to map all
!reduced coordinates xred back into [0,1).
!
!Loop on shells in r-space as was done in g-space
 nr=0
 do
   nr=nr+1
   newr=0
!  
   do ir3=-nr,nr
     do ir2=-nr,nr
       do ir1=-nr,nr
         if( abs(ir3)==nr .or. abs(ir2)==nr .or. abs(ir1)==nr&
&         .or. nr==1 )then

           do ia=1,natom
!            Map reduced coordinate xred(mu,ia) into [0,1)
             fraca1=xred(1,ia)-aint(xred(1,ia))+0.5_dp-sign(0.5_dp,xred(1,ia))
             fraca2=xred(2,ia)-aint(xred(2,ia))+0.5_dp-sign(0.5_dp,xred(2,ia))
             fraca3=xred(3,ia)-aint(xred(3,ia))+0.5_dp-sign(0.5_dp,xred(3,ia))
             drdta1=0.0_dp
             drdta2=0.0_dp
             drdta3=0.0_dp
!            OCL SCALAR ! by MM for Fujitsu
             do ib=1,natom
               fracb1=xred(1,ib)-aint(xred(1,ib))+0.5_dp-sign(0.5_dp,xred(1,ib))
               fracb2=xred(2,ib)-aint(xred(2,ib))+0.5_dp-sign(0.5_dp,xred(2,ib))
               fracb3=xred(3,ib)-aint(xred(3,ib))+0.5_dp-sign(0.5_dp,xred(3,ib))
               r1=dble(ir1)+fracb1-fraca1
               r2=dble(ir2)+fracb2-fraca2
               r3=dble(ir3)+fracb3-fraca3
               rsq=rmet(1,1)*r1*r1+rmet(2,2)*r2*r2+rmet(3,3)*r3*r3+&
&               2.0_dp*(rmet(2,1)*r2*r1+rmet(3,2)*r3*r2+rmet(3,1)*r1*r3)

!              Avoid zero denominators in 'term':
               if (rsq>=1.0d-24) then

!                Note: erfc(8) is about 1.1e-29,
!                so do not bother with larger arg.
!                Also: exp(-64) is about 1.6e-28,
!                so do not bother with larger arg**2 in exp.
                 term=0._dp
                 if (eta*rsq<64.0_dp) then
                   newr=1
                   rmagn=sqrt(rsq)
                   arg=reta*rmagn
!                  derfc is the real(dp) complementary error function
                   !call derfc(derfc_arg,arg)
                   derfc_arg = erfc(arg)
                   term=derfc_arg/rmagn
                   sumr=sumr+zion(typat(ia))*zion(typat(ib))*term
                   term=zion(typat(ia))*zion(typat(ib))*&
&                   (term+fac*exp(-eta*rsq))/rsq
!                  Length scale grads now handled with stress tensor in ewald2
                   r1a1d=rmet(1,1)*r1+rmet(1,2)*r2+rmet(1,3)*r3
                   r2a2d=rmet(2,1)*r1+rmet(2,2)*r2+rmet(2,3)*r3
                   r3a3d=rmet(3,1)*r1+rmet(3,2)*r2+rmet(3,3)*r3
!                  Compute terms related to coordinate gradients
                   drdta1=drdta1+term*r1a1d
                   drdta2=drdta2+term*r2a2d
                   drdta3=drdta3+term*r3a3d
                 end if

!                End avoid zero denominators in'term'
               end if

!              end loop over ib:
             end do

             grewtn(1,ia)=grewtn(1,ia)+drdta1
             grewtn(2,ia)=grewtn(2,ia)+drdta2
             grewtn(3,ia)=grewtn(3,ia)+drdta3

!            end loop over ia:
           end do

!          end triple loop over real space points and associated condition of new shell
         end if
       end do
     end do
   end do

!  Check if new shell must be calculated
   if(newr==0) exit

!  End loop on nr (new shells). Note that there is an exit within the loop
 end do
!
 sumr=0.5_dp*sumr
 fac=pi*ch**2/(2.0_dp*eta*ucvol)

!Finally assemble Ewald energy, eew
 eew=sumg+sumr-chsq*reta/sqrt(pi)-fac

!DEBUG
!write(std_out,*)'eew=sumg+sumr-chsq*reta/sqrt(pi)-fac'
!write(std_out,*)eew,sumg,sumr,chsq*reta/sqrt(pi),fac
!ENDDEBUG

!Length scale grads handled with stress tensor, ewald2

!Output the final values of ng and nr
! write(message, '(a,a,i4,a,i4)' )ch10,&
!& ' ewald : nr and ng are ',nr,' and ',ng
! call wrtout(std_out,message,'COLL')
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
!! PARENTS
!!      stress
!!
!! CHILDREN
!!      derfc,matr3inv
!!
!! SOURCE
!Arguments ------------------------------------
!scalars
 real(dp), parameter :: two_pi = 2*pi
 integer,intent(in) :: natom,ntypat
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: gmet(3,3),rmet(3,3),rprimd(3,3),gprimd(3,3),xred(3,natom)
 real(dp),intent(in) :: zion(ntypat)
 real(dp),intent(out) :: stress(6)

!Local variables-------------------------------
!scalars
 integer :: ia,ib,ig1,ig2,ig3,ir1,ir2,ir3,newg,newr,ng,nr
 real(dp) :: arg1,arg2,arg3,ch,dderfc,derfc_arg,direct,eta,fac,fraca1
 real(dp) :: fraca2,fraca3,fracb1,fracb2,fracb3,g1,g2,g3,gsq,r1,r1c,r2,r2c
 real(dp) :: r3,r3c,recip,reta,rmagn,rsq,summi,summr,t1,t2,t3,t4,t5,t6,term1
 real(dp) :: term2,term3,term4
!arrays
 real(dp) :: strg(6),strr(6)

! *************************************************************************

!Add up total charge and sum of charge^2 in cell
 ch=0._dp
 do ia=1,natom
   ch=ch+zion(typat(ia))
 end do

!Compute eta, the Ewald summation convergence parameter,
!for approximately optimized summations:
 direct=rmet(1,1)+rmet(1,2)+rmet(1,3)+rmet(2,1)+&
& rmet(2,2)+rmet(2,3)+rmet(3,1)+rmet(3,2)+rmet(3,3)
 recip=gmet(1,1)+gmet(1,2)+gmet(1,3)+gmet(2,1)+&
& gmet(2,2)+gmet(2,3)+gmet(3,1)+gmet(3,2)+gmet(3,3)
!Here, a bias is introduced, because G-space summation scales
!better than r space summation !
 eta=pi*200.0_dp/33.0_dp*sqrt(1.69_dp*recip/direct)

 fac=pi**2/eta

!Conduct reciprocal space summations
 strg(1:6)=0.0_dp

!Sum over G space, done shell after shell until all
!contributions are too small
 ng=0
 do
   ng=ng+1
   newg=0

   do ig3=-ng,ng
     do ig2=-ng,ng
       do ig1=-ng,ng

!        Exclude shells previously summed over
         if(abs(ig1)==ng .or. abs(ig2)==ng .or. abs(ig3)==ng&
&         .or. ng==1 ) then

!          Compute Cartesian components of each G
           g1=gprimd(1,1)*ig1+gprimd(1,2)*ig2+gprimd(1,3)*ig3
           g2=gprimd(2,1)*ig1+gprimd(2,2)*ig2+gprimd(2,3)*ig3
           g3=gprimd(3,1)*ig1+gprimd(3,2)*ig2+gprimd(3,3)*ig3
!          Compute |G|^2 (no pi factors)
           gsq=(g1**2+g2**2+g3**2)

!          skip g=0:
           if (gsq>1.0d-20) then
             arg1=fac*gsq

!            larger arg1 gives 0 contribution because of exp(-arg1)
             if (arg1<=80._dp) then
!              When any term contributes then include next shell
               newg=1
               term1=exp(-arg1)/arg1
               summr = 0.0_dp
               summi = 0.0_dp
               do ia=1,natom
                 arg2=two_pi*(ig1*xred(1,ia)+ig2*xred(2,ia)+ig3*xred(3,ia))
!                Sum real and imaginary parts (avoid complex variables)
                 summr=summr+zion(typat(ia))*cos(arg2)
                 summi=summi+zion(typat(ia))*sin(arg2)
               end do

!              Avoid underflow error messages
               if (abs(summr)<1.d-16) summr=0.0_dp
               if (abs(summi)<1.d-16) summi=0.0_dp

               term2=(2._dp/gsq)*(1._dp+arg1)
               t1=term2*g1*g1-1._dp
               t2=term2*g2*g2-1._dp
               t3=term2*g3*g3-1._dp
               t4=term2*g2*g3
               t5=term2*g1*g3
               t6=term2*g1*g2
               term3=term1*(summr*summr+summi*summi)
               strg(1)=strg(1)+t1*term3
               strg(2)=strg(2)+t2*term3
               strg(3)=strg(3)+t3*term3
               strg(4)=strg(4)+t4*term3
               strg(5)=strg(5)+t5*term3
               strg(6)=strg(6)+t6*term3

!              End condition not being larger than 80.0
             end if

!            End skip g=0
           end if

!          End triple loop and condition of new shell
         end if
       end do
     end do
   end do

!  Check if new shell must be calculated
   if (newg==0) exit

!  End loop on new shell. Note that there is an "exit" instruction within the loop
 end do


!Conduct real space summations
 reta=sqrt(eta)
 strr(1:6)=0.0_dp

!Loop on shells in r-space as was done in g-space
 nr=0
 do
   nr=nr+1
   newr=0

   do ir3=-nr,nr
     do ir2=-nr,nr
       do ir1=-nr,nr
         if( abs(ir3)==nr .or. abs(ir2)==nr .or. abs(ir1)==nr&
&         .or. nr==1 )then

           do ia=1,natom
!            Convert reduced atomic coordinates to [0,1)
             fraca1=xred(1,ia)-aint(xred(1,ia))+0.5_dp-sign(0.5_dp,xred(1,ia))
             fraca2=xred(2,ia)-aint(xred(2,ia))+0.5_dp-sign(0.5_dp,xred(2,ia))
             fraca3=xred(3,ia)-aint(xred(3,ia))+0.5_dp-sign(0.5_dp,xred(3,ia))
             do ib=1,natom
               fracb1=xred(1,ib)-aint(xred(1,ib))+0.5_dp-sign(0.5_dp,xred(1,ib))
               fracb2=xred(2,ib)-aint(xred(2,ib))+0.5_dp-sign(0.5_dp,xred(2,ib))
               fracb3=xred(3,ib)-aint(xred(3,ib))+0.5_dp-sign(0.5_dp,xred(3,ib))
               r1=ir1+fracb1-fraca1
               r2=ir2+fracb2-fraca2
               r3=ir3+fracb3-fraca3
!              Convert from reduced to cartesian coordinates
               r1c=rprimd(1,1)*r1+rprimd(1,2)*r2+rprimd(1,3)*r3
               r2c=rprimd(2,1)*r1+rprimd(2,2)*r2+rprimd(2,3)*r3
               r3c=rprimd(3,1)*r1+rprimd(3,2)*r2+rprimd(3,3)*r3
!              Compute |r|^2
               rsq=r1c**2+r2c**2+r3c**2
               rmagn=sqrt(rsq)

!              Avoid zero denominators in 'term':
               if (rmagn>=1.0d-12) then

!                Note: erfc(8) is about 1.1e-29,
!                so do not bother with larger arg.
!                Also: exp(-64) is about 1.6e-28,
!                so do not bother with larger arg**2 in exp.
                 arg3=reta*rmagn
                 if (arg3<8.0_dp) then
                   newr=1
!                  derfc computes the complementary error function
!                  dderfc is the derivative of the complementary error function
                   dderfc=(-2/sqrt(pi))*exp(-eta*rsq)
                   !call derfc(derfc_arg,arg3)
                   derfc_arg = erfc(arg3)
                   term3=dderfc-derfc_arg/arg3
                   term4=zion(typat(ia))*zion(typat(ib))*term3
                   strr(1)=strr(1)+term4*r1c*r1c/rsq
                   strr(2)=strr(2)+term4*r2c*r2c/rsq
                   strr(3)=strr(3)+term4*r3c*r3c/rsq
                   strr(4)=strr(4)+term4*r2c*r3c/rsq
                   strr(5)=strr(5)+term4*r1c*r3c/rsq
                   strr(6)=strr(6)+term4*r1c*r2c/rsq
!                  End the condition of not being to large
                 end if

!                End avoid zero denominator
               end if

!              End loop over ib:
             end do

!            End loop over ia:
           end do

!          End triple loop overs real space points, and associated new shell condition
         end if
       end do
     end do
   end do

!  Check if new shell must be calculated
   if(newr==0) exit

!  End loop on new shells
 end do

!Finally assemble stress tensor coming from Ewald energy, stress
!(note division by unit cell volume in accordance with definition
!found in Nielsen and Martin, Phys. Rev. B 32, 3792 (1985).)

 fac = pi/(2._dp*ucvol*eta)
 stress(1)=(0.5_dp*reta*strr(1)+fac*(strg(1)+(ch**2)))/ucvol
 stress(2)=(0.5_dp*reta*strr(2)+fac*(strg(2)+(ch**2)))/ucvol
 stress(3)=(0.5_dp*reta*strr(3)+fac*(strg(3)+(ch**2)))/ucvol
 stress(4)=(0.5_dp*reta*strr(4)+fac*strg(4))/ucvol
 stress(5)=(0.5_dp*reta*strr(5)+fac*strg(5))/ucvol
 stress(6)=(0.5_dp*reta*strr(6)+fac*strg(6))/ucvol

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
! Calculates the Coulombic energy E as a direc sum. It is very slow, but simple
! to implement, so it is used for testing correctness of more advanced methods.
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
forces = -forces
stress = -stress * ucvol
end subroutine

end module
