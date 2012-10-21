!=========================================================================================
!Copyright (c) 2009, The Regents of the University of Massachusetts, Amherst.
!Developed by E. Polizzi
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without modification, 
!are permitted provided that the following conditions are met:
!
!1. Redistributions of source code must retain the above copyright notice, this list of conditions 
!   and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
!   and the following disclaimer in the documentation and/or other materials provided with the distribution.
!3. Neither the name of the University nor the names of its contributors may be used to endorse or promote
!    products derived from this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
!BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
!ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
!EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
!LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
!IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!==========================================================================================





subroutine sset_point_gauss_legendre(nbe,e,xe,we)
  !  Purpose 
  !  =======
  !  GAUSS-LEGENDRE Quadrature ROUTINE- Return weight and coordinate of the e^th node from nbe total # nodes
  ! 
  !  SINGLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  nbe        (input)        INTEGER: total # nodes for Gauss quadrature
  !  e          (input)        INTEGER: e^th node
  !  xe         (output)       REAL SINGLE PRECISION:  Gauss coordinate for node e 
  !  we         (output)       REAL SINGLE PRECISION:  Gauss weight for node e
  !=====================================================================
  ! Eric Polizzi 2009
  ! ====================================================================
  implicit none
  integer :: nbe,e
  real :: xe,we
!!!!!!!!!!
  double precision ::dxe,dwe

  dxe=dble(xe)
  dwe=dble(we)

  call  dset_point_gauss_legendre(nbe,e,dxe,dwe)

  xe=real(dxe)
  we=real(dwe)

end subroutine sset_point_gauss_legendre



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine dset_point_gauss_legendre(nbe,e,xe,we)
  !  Purpose 
  !  =======
  !  GAUSS-LEGENDRE Quadrature ROUTINE- Return weight and coordinate of the e^th node from nbe total # nodes
  ! 
  !  DOUBLE PRECISION version  
  !
  !  Arguments
  !  =========
  !
  !  nbe        (input)        INTEGER: total # nodes for Gauss quadrature
  !  e          (input)        INTEGER: e^th node
  !  xe         (output)       REAL DOUBLE PRECISION:  Gauss coordinate for node e 
  !  we         (output)       REAL DOUBLE PRECISION:  Gauss weight for node e
  !=====================================================================
  ! Eric Polizzi 2009
  ! ====================================================================
  implicit none
  integer :: nbe,e
  double precision :: xe,we
!!!!!!!!!!

  !! nbe=3
  select case(nbe)

  case(3)

     select case(e)
     case(1)
        we=8.0d0/9.0d0
        xe=0.0d0
     case(2)
        we=5.0d0/9.0d0
        xe=sqrt(3.0d0/5.0d0)
     case(3)
        we=5.0d0/9.0d0 
        xe=-sqrt(3.0d0/5.0d0)
     end select

  case(4)

     select case(e)

     case(1,2)
        we=(18.0d0+sqrt(30.0d0))/36.0d0
        !6.52145154862546142644d-01
        xe=sqrt((3.0d0-2.0d0*sqrt(6.0d0/5.0d0))/7.0d0)
        !3.39981043584856264792d-01
     case(3,4)
        we=(18.0d0-sqrt(30.0d0))/36.0d0
        ! 3.47854845137453857383d-01
        xe=sqrt((3.0d0+2.0d0*sqrt(6.0d0/5.0d0))/7.0d0)
        !8.61136311594052575248d-01

     end select

     if (ibclr(e,0)==e) xe=-xe


  case(5)

     select case(e)
     case(1)
        we=128.0d0/225.0d0
        xe=0.0d0
     case(2,3)
        !      we=9.06179845938663992811d-01
        we=(322.0d0+13.0d0*sqrt(70.0d0))/900.0d0
        !     xe= 5.38469310105683091018d-01 
        xe=(1.0d0/3.0d0)*sqrt(5.0d0-2.0d0*sqrt(10.0d0/7.0d0))
        if (e==3) xe=-xe
     case(4,5)
        !     we=2.36926885056189087515d-01
        we=(322.0d0-13.0d0*sqrt(70.0d0))/900.0d0
        !        xe=9.06179845938663992811d-01
        xe=(1.0d0/3.0d0)*sqrt(5.0d0+2.0d0*sqrt(10.0d0/7.0d0))
        if (e==5) xe=-xe
     end select





  case(6)


     select case(e)
     case(1,2)
        xe=6.61209386466264513688d-01
        we= 3.60761573048138607569d-01
     case(3,4)
        xe=2.38619186083196908630d-01
        we=4.67913934572691047389d-01
     case(5,6)
        xe=9.32469514203152027832d-01
        we=1.71324492379170345043d-01
     end select

     if (ibclr(e,0)==e) xe=-xe




  case(8)


     select case(e)
     case(1,2)
        xe=1.83434642495649804936d-01
        we=3.62683783378361982976d-01
     case(3,4)
        xe=5.25532409916328985830d-01
        we=3.13706645877887287338d-01
     case(5,6)
        xe=7.96666477413626739567d-01
        we=2.22381034453374470546d-01
     case(7,8)
        xe=9.60289856497536231661d-01
        we=1.01228536290376259154d-01
     end select

     if (ibclr(e,0)==e) xe=-xe

  case(10)

     select case(e)
     case(1,2)
        we=2.95524224714752870187d-01
        xe=1.48874338981631210881d-01
     case(3,4)
        we=2.69266719309996355105d-01
        xe=4.33395394129247190794d-01
     case(5,6)
        we=2.19086362515982044000d-01
        xe=6.79409568299024406207d-01
     case(7,8)
        we=1.49451349150580593150d-01
        xe=8.65063366688984510759d-01
     case(9,10)
        we=6.66713443086881375920d-02
        xe=9.73906528517171720066d-01
     end select

     if (ibclr(e,0)==e) xe=-xe


  case(12)

     select case(e)
     case(1,2)
        we=2.49147045813402785006d-01
        xe=1.25233408511468915478d-01
     case(3,4)
        we=2.33492536538354808758d-01
        xe=3.67831498998180193757d-01
     case(5,6)
        we=2.03167426723065921743d-01
        xe=5.87317954286617447312d-01
     case(7,8)
        we=1.60078328543346226338d-01
        xe=7.69902674194304687059d-01
     case(9,10)
        we=1.06939325995318430960d-01
        xe=9.04117256370474856682d-01
     case(11,12)
        we=4.71753363865118271952d-02
        xe=9.81560634246719250712d-01


     end select

     if (ibclr(e,0)==e) xe=-xe




  case(16)


     select case(e)
     case(1,2)
        we=1.89450610455068496287d-01
        xe=9.50125098376374401877d-02
     case(3,4)
        we=1.82603415044923588872d-01
        xe=2.81603550779258913231d-01
     case(5,6)
        we=1.69156519395002538183d-01
        xe=4.58016777657227386350d-01
     case(7,8)
        we=1.49595988816576732080d-01
        xe=6.17876244402643748452d-01
     case(9,10)
        we=1.24628971255533872056d-01
        xe=7.55404408355003033891d-01
     case(11,12)
        we=9.51585116824927848073d-02
        xe=8.65631202387831743866d-01
     case(13,14)
        we=6.22535239386478928628d-02
        xe=9.44575023073232576090d-01
     case(15,16)
        we=2.71524594117540948514d-02
        xe=9.89400934991649932601d-01
     end select

     if (ibclr(e,0)==e) xe=-xe



  case(20)


     select case(e)
     case(1,2)
        we=1.52753387130725850699d-01
        xe=7.65265211334973337513d-02
     case(3,4)
        we=1.49172986472603746785d-01
        xe=2.27785851141645078076d-01
     case(5,6)
        we=1.42096109318382051326d-01
        xe=3.73706088715419560662d-01
     case(7,8)
        we=1.31688638449176626902d-01
        xe=5.10867001950827097985d-01
     case(9,10)
        we=1.18194531961518417310d-01
        xe=6.36053680726515025467d-01  
     case(11,12)
        we=1.01930119817240435039d-01
        xe=7.46331906460150792634d-01
     case(13,14)
        we=8.32767415767047487264d-02
        xe=8.39116971822218823420d-01 
     case(15,16)
        we=6.26720483341090635663d-02
        xe=9.12234428251325905857d-01
     case(17,18)
        we=4.06014298003869413320d-02 
        xe=9.63971927277913791287d-01
     case(19,20)
        we=1.76140071391521183115d-02
        xe=9.93128599185094924776d-01

     end select

     if (ibclr(e,0)==e) xe=-xe




  case(24)


     select case(e)
     case(1,2)
        we=1.27938195346752156976d-01
        xe=6.40568928626056260827d-02
     case(3,4)
        we=1.25837456346828296117d-01
        xe=1.91118867473616309153d-01
     case(5,6)
        we=1.21670472927803391202d-01
        xe=3.15042679696163374398d-01
     case(7,8)
        we=1.15505668053725601353d-01
        xe=4.33793507626045138478d-01
     case(9,10)
        we=1.07444270115965634785d-01
        xe=5.45421471388839535649d-01
     case(11,12)
        we=9.76186521041138882720d-02
        xe=6.48093651936975569268d-01
     case(13,14)
        we=8.61901615319532759152d-02
        xe=7.40124191578554364260d-01
     case(15,16)
        we=7.33464814110803057346d-02
        xe=8.20001985973902921981d-01
     case(17,18)
        we=5.92985849154367807461d-02
        xe=8.86415527004401034190d-01
     case(19,20)
        we=4.42774388174198061695d-02
        xe=9.38274552002732758539d-01
     case(21,22)
        we=2.85313886289336631809d-02
        xe=9.74728555971309498199d-01
     case(23,24)
        we=1.23412297999871995469d-02
        xe=9.95187219997021360195d-01
     end select

     if (ibclr(e,0)==e) xe=-xe


  case(32)


     select case(e)
     case(1,2)
        we=9.65400885147278005666d-02
        xe=4.83076656877383162364d-02
     case(3,4)
        we=9.56387200792748594185d-02
        xe=1.44471961582796493484d-01
     case(5,6)
        we=9.38443990808045656367d-02
        xe=2.39287362252137074544d-01
     case(7,8)
        we=9.11738786957638847129d-02
        xe=3.31868602282127649782d-01
     case(9,10)
        we=8.76520930044038111450d-02
        xe=4.21351276130635345353d-01
     case(11,12)
        we=8.33119242269467552223d-02
        xe=5.06899908932229390044d-01
     case(13,14)
        we=7.81938957870703064685d-02
        xe=5.87715757240762329066d-01
     case(15,16)
        we=7.23457941088485062287d-02
        xe=6.63044266930215200960d-01
     case(17,18)
        we=6.58222227763618468406d-02
        xe=7.32182118740289680412d-01
     case(19,20)
        we=5.86840934785355471448d-02
        xe=7.94483795967942406965d-01
     case(21,22)
        we=5.09980592623761761959d-02
        xe=8.49367613732569970160d-01
     case(23,24)
        we=4.28358980222266806557d-02
        xe=8.96321155766052123971d-01
     case(25,26)
        we=3.42738629130214331033d-02
        xe=9.34906075937739689159d-01
     case(27,28)
        we=2.53920653092620594561d-02
        xe=9.64762255587506430761d-01
     case(29,30)
        we=1.62743947309056706058d-02
        xe=9.85611511545268335400d-01
     case(31,32)
        we=7.01861000947009660028d-03
        xe=9.97263861849481563534d-01
     end select

     if (ibclr(e,0)==e) xe=-xe





  case(40)

     select case(e)
     case(1,2)
        we=7.75059479784248112668d-02
        xe=3.87724175060508219329d-02
     case(3,4)
        we=7.70398181642479655914d-02
        xe=1.16084070675255208481d-01
     case(5,6)
        we=7.61103619006262423723d-02
        xe=1.92697580701371099719d-01
     case(7,8)
        we=7.47231690579682641980d-02
        xe=2.68152185007253681152d-01
     case(9,10)
        we=7.28865823958040590609d-02
        xe=3.41994090825758473008d-01
     case(11,12)
        we=7.06116473912867796979d-02
        xe=4.13779204371605001525d-01
     case(13,14)
        we=6.79120458152339038265d-02
        xe=4.83075801686178712903d-01
     case(15,16)
        we=6.48040134566010380719d-02
        xe=5.49467125095128202056d-01
     case(17,18)
        we=6.13062424929289391679d-02    
        xe=6.12553889667980237972d-01
     case(19,20)
        we=5.74397690993915513665d-02
        xe=6.71956684614179548364d-01
     case(21,22)
        we=5.32278469839368243566d-02
        xe=7.27318255189927103277d-01
     case(23,24)
        we=4.86958076350722320604d-02
        xe=7.78305651426519387712d-01
     case(25,26)
        we=4.38709081856732719923d-02
        xe=8.24612230833311663197d-01
     case(27,28)
        we=3.87821679744720176413d-02
        xe=8.65959503212259503824d-01
     case(29,30)
        we=3.34601952825478473933d-02
        xe=9.02098806968874296732d-01
     case(31,32)
        we=2.79370069800234010984d-02
        xe=9.32812808278676533383d-01
     case(33,34)
        we=2.22458491941669572615d-02
        xe=9.57916819213791655824d-01
     case(35,36)
        we=1.64210583819078887131d-02
        xe=9.77259949983774262679d-01
     case(37,38)
        xe=9.90726238699457006464d-01
        we=1.04982845311528136146d-02
     case(39,40)
        we=4.52127709853319125846d-03
        xe=9.98237709710559200369d-01
     end select

     if (ibclr(e,0)==e) xe=-xe



  case(48)


     select case(e)
     case(1,2)
        we=6.47376968126839225006d-02
        xe=3.23801709628693620343d-02
     case(3,4)
        we=6.44661644359500822082d-02
        xe=9.70046992094626989322d-02
     case(5,6)
        we=6.39242385846481866207d-02
        xe=1.61222356068891718055d-01
     case(7,8)
        we=6.31141922862540256548d-02
        xe=2.24763790394689061224d-01
     case(9,10)
        we=6.20394231598926639029d-02
        xe=2.87362487355455576728d-01
     case(11,12)
        we=6.07044391658938800517d-02
        xe=3.48755886292160738148d-01
     case(13,14)
        we=5.91148396983956357477d-02
        xe=4.08686481990716729925d-01
     case(15,16)
        we=5.72772921004032157044d-02
        xe=4.66902904750958404535d-01
     case(17,18)
        we=5.51995036999841628676d-02
        xe=5.23160974722233033658d-01
     case(19,20)
        we=5.28901894851936670964d-02
        xe=5.77224726083972703838d-01
     case(21,22)
        we=5.03590355538544749590d-02
        xe=6.28867396776513624013d-01
     case(23,24)
        we=4.76166584924904748267d-02
        xe=6.77872379632663905208d-01
     case(25,26)
        we=4.46745608566942804201d-02
        xe=7.24034130923814654658d-01
     case(27,28)
        we=4.15450829434647492133d-02
        xe=7.67159032515740339276d-01
     case(29,30)
        we=3.82413510658307063158d-02
        xe=8.07066204029442627087d-01
     case(31,32)
        we=3.47772225647704388909d-02
        xe=8.43588261624393530704d-01
     case(33,34)
        we=3.11672278327980889025d-02
        xe=8.76572020274247885885d-01
     case(35,36)
        we=2.74265097083569482001d-02
        xe=9.05879136715569672805d-01
     case(37,38)
        we=2.35707608393243791410d-02
        xe=9.31386690706554333107d-01
     case(39,40)
        we=1.96161604573555278142d-02
        xe=9.52987703160430860724d-01
     case(41,42)
        we=1.55793157229438487279d-02
        xe=9.70591592546247250472d-01
     case(43,44)
        we=1.14772345792345394895d-02
        xe=9.84124583722826857765d-01
     case(45,46)
        we=7.32755390127626210220d-03
        xe=9.93530172266350757526d-01
     case(47,48)
        we=3.15334605230583863260d-03
        xe=9.98771007252426118580d-01

     end select

     if (ibclr(e,0)==e) xe=-xe



  end select!!!!!!!!!!!!!!!!!!!!!!! nbe


end subroutine dset_point_gauss_legendre
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



