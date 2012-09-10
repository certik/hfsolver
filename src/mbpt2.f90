module mbpt2
use types, only: dp
use scf, only: ijkl2intindex
implicit none

contains

real(dp) function term1(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, r, s
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do r = Noccupied+1, size(lam)
            do s = Noccupied+1, size(lam)
                E2 = E2 + (  &
                        + (-1)**(2+2) * 2**2 / 2._dp  &
                      * moint2(ijkl2intindex(a, r, b, s))  &
                      * moint2(ijkl2intindex(r, a, s, b))  &
                        + (-1)**(1+2) * 2**1 / 2._dp  &
                      * moint2(ijkl2intindex(a, r, b, s))  &
                      * moint2(ijkl2intindex(s, a, r, b))  &
                    ) / (  &
                    (-lam(r)-lam(s)+lam(a)+lam(b))  &
                    )
            end do
        end do
    end do
end do
end function

end module
