program test_bessel
use types, only: dp
use constants, only: pi
use special_functions, only: Fm
use utils, only: assert, init_random
implicit none
integer :: i, n
real(dp) :: r
real(dp), allocatable :: yf(:), yr(:), x(:)
real(dp) :: t1, t2, t3

call init_random()

n = 10000000
allocate(x(n), yf(n), yr(n))
do i = 1, n
    call random_number(r)
    x(i) = r*10 + 10
end do
call cpu_time(t1)
do i = 1, n
    yf(i) = f(x(i))
end do
call cpu_time(t2)
do i = 1, n
    yr(i) = r10_20(x(i))
end do
call cpu_time(t3)
print *, "abs:", maxval(abs(yf-yr))
print *, "rel:", maxval(abs(yf-yr) / max(abs(yf), abs(yr)))
print *, "time f(r):", t2-t1
print *, "time r(r):", t3-t2
print *, "speedup:", (t2-t1) / (t3-t2)

contains

real(dp) function f(x) result(r)
real(dp) :: x
r = ((105/x**4 + 45/x**2 + 1)*sinh(x) - (105/x**3 + 10/x)*cosh(x)) / exp(x)
end function

real(dp) function r4_10(x) result(r)
real(dp) :: x
        r = (0.000395502959013236968661582656143_dp - &
                0.001434648369704841686633794071_dp*x + &
                0.00248783474583503473135143644434_dp*x**2 - &
                0.00274477921388295929464613063609_dp*x**3 + &
                0.00216275018107657273725589740499_dp*x**4 - &
                0.000236779926184242197820134964535_dp*x**5 + &
                0.0000882030507076791807159699814428_dp*x**6 - &
                4.62078105288798755556136693122e-6_dp*x**7 + &
                8.23671374777791529292655504214e-7_dp*x**8) /&
            (1 + 0.504839286873735708062045336271_dp*x + &
                0.176683950009401712892997268723_dp*x**2 + &
                0.0438594911840609324095487447279_dp*x**3 + &
                0.00829753062428409331123592322788_dp*x**4 + &
                0.00111693697900468156881720995034_dp*x**5 + &
                0.000174719963536517752971223459247_dp*x**6 + &
                7.22885338737473776714257581233e-6_dp*x**7 + &
                1.64737453771748367647332279826e-6_dp*x**8)
end function

real(dp) function r10_20(x) result(r)
real(dp) :: x
        r = (1.49435717183021678294278540018_dp + &
            x*(-1.9954827594990599398954087063_dp + &
            x*(1.19185825369343226912112655137_dp + &
            x*(-0.40866680980235804096143699423_dp + &
            x*(0.0852839860059780325406440673318_dp + &
            (-0.00980617919194154929317057489645_dp + &
            0.000550291361244287676343295379476_dp*x)*x)))))/&
        (1 + x*(0.420439518058743727857466136746_dp + &
            x*(0.144024726914933127664739439568_dp + &
            x*(0.035261250406130055921113600336_dp + &
            x*(0.0349770458351085078647522073879_dp + &
            (-0.00860653991097136433951965579037_dp + &
            0.00110058277850687516223459976889_dp*x)*x)))))
end function

end program
