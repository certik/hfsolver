module pulay
!
!  purpose:  
!       calculate mixing for accelerated convergence in an interative 
!       procedure such as charge-density mixing in scf calculations
!
!  usage:
!       call pulayMix(cdin,cdout,npulay,linearmixing,outputmixing,distance)
!           where
!             cdin()         input cd for the current iteration on input 
!                            and is overwritten by the input cd for the 
!                            next iteration on output
!             cdout()        is the output cd for the current iteration
!             npulay         number of iterations used in cd update  
!             linearmixing   is the linear mixing between input and output
!                            charge densities for the first iteration.  This
!                            is a value between zero and 1, so a value of 
!                            0.1 means mix 10% of output with 90% of input.
!              outputmixing  parameter between zero and 1.  
!                            If present and non-zero, mix the optimal output 
!                            cd with the optimal input cd to form new input.
!                            A value of 0.1 means mix 10% of optimal output
!                            with 90% of optimal input. 
!              distance      output parameter giving the pulay distance 
!
!              
!                
!      The pulay module makes use of two scratch files in order to minimize 
!      internal memory requirements.  
!      These files have hardwired unit numbers: 31,32
!                
!         unit 31:     formatted file containing the iteration number, the 
!                      number of cd files to use in the mixing, and the 
!                      matrix of difference dot-products.
!
!         unit 32:     binary file formatted for direct access with 
!                      record length set to the size of the charge 
!                      density file 
!                      contains charge density and residuals for the 
!                      most recent iterations
!  revised 6/2010, J.E. Pask, LLNL: fixed diffarray bounds in lines 235--241
!  3/2011: removed extraneous "use matrixVector" --J.E.P, LLNL
use types
use lapack, only: dsysv
use utils, only: stop_error

implicit none
private
public pulayMix,linearMix

CONTAINS

subroutine linearMix(cdin,cdout,linearmixing)
real (dp), intent(inout) :: cdin(:)         ! input cd to current iteration
                                            ! overwritten by next input
real (dp), intent(in)    :: cdout(:)        ! output cd from current iteration
real (dp), intent(in)    :: linearmixing    ! mixing parameter - 0.1 means
                                            ! mix 10% or output with input
!
!  check that linear mixing parameter is between 0 and 1 and that input and 
!  output charge densities have the same dimension
!
if(linearmixing < 0._dp .OR. linearmixing > 1._dp) then
   write(6,*) ' error in linearMix: linearmixing out of range'
   write(6,*) 'linearmixing=',linearmixing
   call stop_error('linearMix:linearmixing')
endif

if (size(cdin) .ne. size(cdout) ) then
   write(6,*) ' error in linearMix: mixmatch in cd array dimensions'
   write(6,*) 'size(cdin),size(cdout)=',size(cdin),size(cdout)
   call stop_error('linearMix:cdinout mismatch')
endif

! now just mix them and we're done!   
cdin = cdin + linearmixing*(cdout-cdin)

end subroutine linearMix


subroutine pulayMix(ncd,cdin,cdout,npulay,linearmixing,outputmixing,distance)

integer, intent(in):: ncd                ! dimension of cd arrays
real (dp), intent(inout) :: cdin(ncd)    ! input cd for the current iteration 
                                         ! overwritten by the input cd for the 
                                         ! next iteration
real (dp), intent(in)    :: cdout(ncd)   ! output cd for the current iteration
integer, intent(inout)   :: npulay       ! number of iterations used in pulay
real (dp), intent(in)    :: linearmixing ! linear mixing between input & output
                                         ! charge densities for 1st iteration.
                                         ! 0<=linearmixing<=1: 0.1 means mix 
                                         ! 10% of output with 90% of input.
real (dp), intent(in)    :: outputmixing ! optional parameter in range [0,1].  
                                         ! If present and non-zero, mix the 
                                         ! optimal output cd with the optimal 
                                         ! input cd to form new input.
                                         ! A value of 0.1 means mix 10% of 
                                         ! optimal output with 90% of optimal 
                                         ! input. 
real (dp), intent(out)   :: distance     ! distance parameter for convergence
integer                  :: ipulay       ! iteration number 
integer                  :: nfile1,nfile2! unit file numbers for pulay files
real(dp)                 :: diffdot      ! variable storing dot-product
real (dp), allocatable   :: diffarray(:,:) ! array of cd difference dotproducts
real (dp), allocatable   :: diffarrayold(:,:) ! array of cd difference dotproducts
integer                  :: ndim         ! dimension of cds in pulay files
real (dp), allocatable   :: diffdots(:)  ! dot-prods: diffs with current diff 
integer                  :: npmat        ! dimension of diffs array
integer                  :: jpulay       ! index of current cd in diffs array
integer                  :: i            ! do-loop index
real (dp), allocatable   :: cddiff(:)    ! cd differences read from nfile2
integer                  :: nplin        ! dimension of matrix for linear equations
real (dp), allocatable   :: linearmatrix(:,:) ! Pulay matrix A for linear equation AX=B
real (dp), allocatable   :: linearvector(:)   ! Vector B for for linear equation AX=B
real (dp), allocatable   :: work(:)      ! work array for LAPACK routine DSYSV
integer                  :: lwork        ! dimension of work array
integer, allocatable     :: ipiv(:)      ! work array for LAPACK routine DSYSV
integer                  :: info         ! diagnostic output from LAPACK routine DSYSV
logical                  :: doOutputMix  ! true if pulay mixes input and output
integer                  :: nold         ! dimension of pulay array in nfile1
logical                :: okfile  ! after inquire, true if file exists
! logical                :: okunit  ! true if getlun returns a unit number
integer                :: precl   ! length of record in pulay2.dat file

! check that charge density files are consistent
ndim = size(cdin)

if (ndim .ne. size(cdout) ) then
   write(6,*) ' error in pulayMix: mixmatch in cd array dimensions'
   write(6,*) 'size(cdin),size(cdout)=',ndim,size(cdout)
   call stop_error('pulayMix:cdinout mismatch')
endif

! do linear mixing for npulay=1 case
if(npulay.eq.1) then
   call linearMix(cdin,cdout,linearmixing)
   diffdot = dot_product(cdout-cdin,cdout-cdin)! compute difference dot-product
   distance = sqrt(diffdot)
else
! check if pulay scratch file is open
inquire (unit=31,opened=okfile)

if (.not. okfile) then  ! first iteration: set up and open new pulay files

   nfile1 = 31
   open(unit=nfile1,form='formatted',status='scratch')
   nfile2 = 32
   inquire (iolength=precl) cdin(:)   ! needed for direct access open
   open(unit=nfile2,form='unformatted',status='scratch', &
   &   access='direct',recl=precl)

   ipulay = 1    ! initialize iteration counter

   diffdot = dot_product(cdout-cdin,cdout-cdin)! compute difference dot-product
                                               ! write and close nfile1, nfile2
   write (nfile1,*) ipulay,npulay,1,ndim,precl
   write (nfile1,*) diffdot

   write (nfile2,rec=1) cdin(:)
   write (nfile2,rec=2) cdout(:)-cdin(:)

   distance = sqrt(diffdot)
                                  ! use linear mixing to form new input
   call linearMix(cdin,cdout,linearmixing)

else   ! second and higher iterations 

! assign unit numbers and open existing pulay files

   nfile1 = 31

   ! read number of previous iterations ipulay and number retained npulay
   ! the value of npulay overwrites the value input to this routine 
   ! allocate and read diffarray
   ! increment ipulay by 1 to reflect the current iteration
   rewind (nfile1)
   read (nfile1,*) ipulay,npulay,npmat,ndim,precl
   allocate(diffarrayold(npmat,npmat))
   read (nfile1,*) diffarrayold(:,:)
   ipulay = ipulay+1

   nfile2 = 32

! for each iteration after the first, we add or replace a charge density
! record.  For ipulay =< npulay, we add a new charge density record.  
! For ipulay>npulay, we overwrite the earliest one.  We retain at most 
! npulay charge densities. We perform similar gymnastics in diffarray. 
! jpulay is the index of the charge density we are adding or replacing.
! nold is the dimension of the previous pulay difference dot-product array. 
   npmat  = min(ipulay,npulay)
   jpulay = mod(ipulay-1,npulay)+1
   nold   = min(ipulay-1,npulay)
                   ! check cd size consistency between pulay file and cdin
   if (ndim .ne. size(cdin)) then
      write (6,*) &
        & 'error in pulayMix: cd array and pulay file cds are different sizes'
      write(6,*) 'ndim,size(cdin)=',ndim,size(cdin)
      write(6,*) 'stopping, ipulay=',ipulay
      call stop_error('pulayMix:pulay/cdin mismatch')
   endif
  ! allocate dot-product array for cd diffs and 
  ! take dot-products with all previous diffs in nfile2
   allocate(diffdots(npmat))
   allocate(cddiff(ndim))

   do i =1,npmat
      if(i .eq. jpulay) then    ! insert/add the current cd in nfile2
         write(nfile2,rec=2*i-1) cdin(:)
         write(nfile2,rec=2*i) cdout(:)-cdin(:)
         diffdots(i) = dot_product(cdout(:)-cdin(:),cdout(:)-cdin(:))
      else                      ! read cd difference and form dot_products 
         read(nfile2,rec=2*i) cddiff(:)
         diffdots(i) = dot_product(cdout(:)-cdin(:),cddiff(:))
      endif
   end do
                   ! form new diffarray matrix
   allocate(diffarray(npmat,npmat))
   diffarray(1:nold,1:nold) = diffarrayold(:,:)
   diffarray(jpulay,:) = diffdots(:)
   diffarray(:,jpulay) = diffdots(:)
                   ! write and close nfile1
   rewind (nfile1)
   write (nfile1,*) ipulay, npulay,npmat,ndim,precl
   write (nfile1,*) diffarray(:,:)


! set up linear equation system with constraint on coefficient sum:
!   substitute in for x_n using sum(x_i)=1
!   then we get a set of linear equations to solve for the 
!   coefficients x_1.....x_n-1 given by A(i,j)X(j)=B(i)
!   where A(i,j) = D(i,j) - D(i,n) - D(n,j) + D(n,n)
!   and   B(i)   = D(n,n) - D(n,i)
!   where D is the diffarray matrix

   nplin=npmat-1
   allocate(linearmatrix(nplin,nplin))
   linearmatrix(:,:) = diffarray(1:nplin,1:nplin) + diffarray(npmat,npmat)
   do i =1,nplin
      linearmatrix(:,i) = linearmatrix(:,i) - diffarray(1:nplin,npmat)
      linearmatrix(i,:) = linearmatrix(i,:) - diffarray(npmat,1:nplin)
   end do
   allocate(linearvector(npmat))
   linearvector(1:nplin) = diffarray(npmat,npmat) - diffarray(npmat,1:nplin)
!   write(6,*) 'nplin=',nplin
!   write(6,*) 'linearmatrix=',linearmatrix
!   write(6,*) 'linearvector=',linearvector(1:nplin)

! solve using LAPACK routine
   allocate(ipiv(nplin))
   lwork = 5*nplin*nplin
   allocate(work(lwork))

   call dsysv('U', nplin, 1, linearmatrix, nplin, ipiv, &
             &               linearvector, nplin, work, lwork, info )
   linearvector(npmat) = 1._dp - sum(linearvector(1:nplin))
!
!  Explicitly evaluate the residual
!
   distance = sqrt(dot_product(linearvector,matmul(diffarray,linearvector)))
  if(info .ne. 0) write(6,*) 'info=',info

!   write(6,*) 'outlinvector=',linearvector
!
! make new optimal input cd
   cdin(:) = 0._dp
   doOutputMix = ( outputmixing .gt. 0._dp .and. outputmixing .le. 1._dp)
   do i = 1,npmat
      read (nfile2,rec=2*i-1) cddiff(:)
      cdin(:) = cdin(:) + linearvector(i)*cddiff(:)
      if(doOutputMix) then   ! mix with new optimal output cd, if necessary
         read (nfile2,rec=2*i) cddiff(:)
      ! Note: cdout here is really the difference cdout-cdin,
      ! so this results in outputmixing*cdout + (1-outputmixing)*cdin
         cdin(:) = cdin(:) + outputmixing*linearvector(i)*cddiff(:)
      endif
   end do

   deallocate(diffarrayold, diffarray, diffdots, cddiff)
   deallocate(linearmatrix, linearvector, ipiv, work)
endif

endif
end subroutine pulayMix

end module pulay
