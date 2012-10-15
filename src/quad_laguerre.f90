module quad_laguerre
! Gauss-Laguerre quadrature routines
use types, only: dp
use utils, only: stop_error, str
implicit none
private
public gauss_laguerre_pts, gauss_laguerre_wts

contains

function gauss_laguerre_pts(n) result(xi)
! Returns Gauss-Laguerre points for n-point quadrature
integer, intent(in):: n ! number of Gauss points
real(dp) xi(n) ! Gauss points for n-point quadrature
select case(n)
    case(6)
        xi = [0.222846604179_dp, 0.118893210173e1_dp, 0.299273632659e1_dp, &
            0.577514356905e1_dp, 0.983746741883e1_dp, 0.159828739802e2_dp]
    case(20)
        xi = [7.05398896919887E-02_dp, 3.72126818001611E-01_dp, &
            9.16582102483273E-01_dp, 1.70730653102834E00_dp, &
            2.74919925530943E00_dp, 4.04892531385089E00_dp, &
            5.61517497086162E00_dp, 7.45901745367106E00_dp, &
            9.59439286958110E00_dp, 1.20388025469643E01_dp, &
            1.48142934426307E01_dp, 1.79488955205194E01_dp, &
            2.14787882402850E01_dp, 2.54517027931869E01_dp, &
            2.99325546317006E01_dp, 3.50134342404790E01_dp, &
            4.08330570567286E01_dp, 4.76199940473465E01_dp, &
            5.58107957500639E01_dp, 6.65244165256157E01_dp]
    case(52)
        xi = [.2753976830875507E-01_dp, .1451324895939190E00_dp, &
            .3568014904470186E00_dp, .6627839659547586E00_dp, &
            .1063361679923641E01_dp, .1558899184297198E01_dp, &
            .2149849112008393E01_dp, .2836755125378294E01_dp, &
            .3620254892936698E01_dp, .4501083489774322E01_dp, &
            .5480077344494273E01_dp, .6558178807399577E01_dp, &
            .7736441409772490E01_dp, .9016035891246450E01_dp, &
            .1039825708441186E02_dp, .1188453176140202E02_dp, &
            .1347642756623108E02_dp, .1517566317961300E02_dp, &
            .1698411989071981E02_dp, .1890385478398055E02_dp, &
            .2093711579007395E02_dp, .2308635890067340E02_dp, &
            .2535426790879905E02_dp, .2774377711415792E02_dp, &
            .3025809753002753E02_dp, .3290074725097065E02_dp, &
            .3567558679693149E02_dp, .3858686044992265E02_dp, &
            .4163924485954081E02_dp, .4483790653382173E02_dp, &
            .4818857028192577E02_dp, .5169760127699361E02_dp, &
            .5537210422243494E02_dp, .5922004422319802E02_dp, &
            .6325039552071346E02_dp, .6747332645345497E02_dp, &
            .7190043217804184E02_dp, .7654503134565677E02_dp, &
            .8142254992297490E02_dp, .8655102610563451E02_dp, &
            .9195178728603965E02_dp, .9765037780519004E02_dp, &
            .1036778632153708E03_dp, .1100727197636936E03_dp, &
            .1168836718458836E03_dp, .1241741438275154E03_dp, &
            .1320296369847717E03_dp, .1405708428393328E03_dp, &
            .1499792501042790E03_dp, .1605542735938268E03_dp, &
            .1728700383236378E03_dp, .1884083008266166E03_dp]
    case default
        call stop_error("gauss_laguerre_pts: n = " // str(n) //" not supported")
end select
end function

function gauss_laguerre_wts(n) result(w)
! Returns Gauss-Laguerre weights for n-point quadrature
integer, intent(in):: n ! number of Gauss points
real(dp) w(n) ! Gauss weights for n-point quadrature
select case(n)
    case(6)
        w = [0.458964673950_dp, 0.417000830772_dp, 0.113373382074_dp, &
            0.103991974531e-1_dp, 0.261017202815e-3_dp, 0.898547906430e-6_dp]
    case(20)
        w = [1.68746801851114E-01_dp, 2.91254362006068E-01_dp, &
            2.66686102867001E-01_dp, 1.66002453269507E-01_dp, &
            7.48260646687924E-02_dp, 2.49644173092832E-02_dp, &
            6.20255084457224E-03_dp, 1.14496238647691E-03_dp, &
            1.55741773027812E-04_dp, 1.54014408652249E-05_dp, &
            1.08648636651798E-06_dp, 5.33012090955671E-08_dp, &
            1.75798117905058E-09_dp, 3.72550240251232E-11_dp, &
            4.76752925157819E-13_dp, 3.37284424336244E-15_dp, &
            1.15501433950040E-17_dp, 1.53952214058234E-20_dp, &
            5.28644272556916E-24_dp, 1.65645661249902E-28_dp]
    case(52)
        w = [.6875910298597036E-01_dp, .1423540973986236E00_dp, &
            .1811266787409769E00_dp, .1820580409718219E00_dp, &
            .1546877197634628E00_dp, .1142633872117745E00_dp, &
            .7442476052638809E-01_dp, .4308762935349540E-01_dp, &
            .2227840111470386E-01_dp, .1031755067884288E-01_dp, &
            .4287330009884616E-02_dp, .1600025983279153E-02_dp, &
            .5364883404246371E-03_dp, .1616067225512215E-03_dp, &
            .4371626901141589E-04_dp, .1061204120019816E-04_dp, &
            .2309393970187362E-05_dp, .4499907198830564E-06_dp, &
            .7839262646345511E-07_dp, .1218906919667385E-07_dp, &
            .1688277042283522E-08_dp, .2078462389876053E-09_dp, &
            .2268826333372677E-10_dp, .2189963751640542E-11_dp, &
            .1863530337037668E-12_dp, .1393307645786146E-13_dp, &
            .9119304426764297E-15_dp, .5203589029695487E-16_dp, &
            .2576930144884706E-17_dp, .1101990742164216E-18_dp, &
            .4046717938672477E-20_dp, .1268159429347318E-21_dp, &
            .3367962147504498E-23_dp, .7521222948802292E-25_dp, &
            .1399938482573673E-26_dp, .2150241263967673E-28_dp, &
            .2694448929466313E-30_dp, .2718730293186383E-32_dp, &
            .2175622513648559E-34_dp, .1356483523043116E-36_dp, &
            .6452804249835994E-39_dp, .2283907890830233E-41_dp, &
            .5833769441130827E-44_dp, .1035516957780775E-46_dp, &
            .1217736979511202E-49_dp, .8915957957460291E-53_dp, &
            .3738606192685472E-56_dp, .7981307307062951E-60_dp, &
            .7272135627891910E-64_dp, .2116837891015250E-68_dp, &
            .1134265318578078E-73_dp, .2750867140614681E-80_dp]
    case default
        call stop_error("gauss_laguerre_wts: n = " // str(n) //" not supported")
end select
end function

end module