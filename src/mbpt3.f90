module mbpt3_mod
use types, only: dp
use scf, only: ijkl2intindex
implicit none

contains

real(dp) function term1(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, r, s, t
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do r = Noccupied+1, size(lam)
                do s = Noccupied+1, size(lam)
                    do t = Noccupied+1, size(lam)
                        E2 = E2 + (  &
                                + (-1)**(3+3) * 2**3  &
                              * moint2(ijkl2intindex(a, r, b, s))  &
                              * moint2(ijkl2intindex(r, a, c, t))  &
                              * moint2(ijkl2intindex(s, b, t, c))  &
                                + (-1)**(2+3) * 2**2  &
                              * moint2(ijkl2intindex(a, r, b, s))  &
                              * moint2(ijkl2intindex(r, a, c, t))  &
                              * moint2(ijkl2intindex(t, b, s, c))  &
                                + (-1)**(2+3) * 2**2  &
                              * moint2(ijkl2intindex(a, r, b, s))  &
                              * moint2(ijkl2intindex(c, a, r, t))  &
                              * moint2(ijkl2intindex(s, b, t, c))  &
                                + (-1)**(1+3) * 2**1  &
                              * moint2(ijkl2intindex(a, r, b, s))  &
                              * moint2(ijkl2intindex(c, a, r, t))  &
                              * moint2(ijkl2intindex(t, b, s, c))  &
                                + (-1)**(2+3) * 2**2  &
                              * moint2(ijkl2intindex(b, r, a, s))  &
                              * moint2(ijkl2intindex(r, a, c, t))  &
                              * moint2(ijkl2intindex(s, b, t, c))  &
                                + (-1)**(1+3) * 2**1  &
                              * moint2(ijkl2intindex(b, r, a, s))  &
                              * moint2(ijkl2intindex(r, a, c, t))  &
                              * moint2(ijkl2intindex(t, b, s, c))  &
                                + (-1)**(1+3) * 2**1  &
                              * moint2(ijkl2intindex(b, r, a, s))  &
                              * moint2(ijkl2intindex(c, a, r, t))  &
                              * moint2(ijkl2intindex(s, b, t, c))  &
                                + (-1)**(2+3) * 2**2  &
                              * moint2(ijkl2intindex(b, r, a, s))  &
                              * moint2(ijkl2intindex(c, a, r, t))  &
                              * moint2(ijkl2intindex(t, b, s, c))  &
                            ) / (  &
                            (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                            (-lam(s)-lam(t)+lam(b)+lam(c))  &
                            )
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term2(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        E2 = E2 + (  &
                                + (-1)**(2+4) * 2**2 / 2._dp  &
                              * moint2(ijkl2intindex(a, r, b, s))  &
                              * moint2(ijkl2intindex(c, a, d, b))  &
                              * moint2(ijkl2intindex(r, c, s, d))  &
                                + (-1)**(1+4) * 2**1 / 2._dp  &
                              * moint2(ijkl2intindex(a, r, b, s))  &
                              * moint2(ijkl2intindex(c, a, d, b))  &
                              * moint2(ijkl2intindex(s, c, r, d))  &
                            ) / (  &
                            (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                            (-lam(r)-lam(s)+lam(c)+lam(d))  &
                            )
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term3(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do r = Noccupied+1, size(lam)
            do s = Noccupied+1, size(lam)
                do t = Noccupied+1, size(lam)
                    do u = Noccupied+1, size(lam)
                        E2 = E2 + (  &
                                + (-1)**(2+2) * 2**2 / 2._dp  &
                              * moint2(ijkl2intindex(a, r, b, s))  &
                              * moint2(ijkl2intindex(r, t, s, u))  &
                              * moint2(ijkl2intindex(t, a, u, b))  &
                                + (-1)**(1+2) * 2**1 / 2._dp  &
                              * moint2(ijkl2intindex(a, r, b, s))  &
                              * moint2(ijkl2intindex(r, t, s, u))  &
                              * moint2(ijkl2intindex(u, a, t, b))  &
                            ) / (  &
                            (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                            (-lam(t)-lam(u)+lam(a)+lam(b))  &
                            )
                    end do
                end do
            end do
        end do
    end do
end do
end function

end module
