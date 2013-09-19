module mbpt4_mod
use types, only: dp
use scf, only: ijkl2intindex
implicit none

contains

real(dp) function term1(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(4+4) * 2**4  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(s, b, d, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(s, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(d, b, s, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(d, b, s, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(s, b, d, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(s, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, s, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, s, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(s, b, d, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(s, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(d, b, s, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(d, b, s, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(s, b, d, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(s, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, s, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, s, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(s)-lam(t)+lam(b)+lam(c)) *  &
                                    (-lam(t)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
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
integer :: a, b, c, d, e, r, s, t
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do e = 1, Noccupied
                    do r = Noccupied+1, size(lam)
                        do s = Noccupied+1, size(lam)
                            do t = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(r, b, e, c))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(r, b, e, c))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(e, b, r, c))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(3+5) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(e, b, r, c))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(d, a, c, t))  &
                                      * moint2(ijkl2intindex(r, b, e, c))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(d, a, c, t))  &
                                      * moint2(ijkl2intindex(r, b, e, c))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(d, a, c, t))  &
                                      * moint2(ijkl2intindex(e, b, r, c))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(d, a, c, t))  &
                                      * moint2(ijkl2intindex(e, b, r, c))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(r, b, e, c))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(3+5) * 2**3  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(r, b, e, c))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(e, b, r, c))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(e, b, r, c))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                        + (-1)**(3+5) * 2**3  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(d, a, c, t))  &
                                      * moint2(ijkl2intindex(r, b, e, c))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(d, a, c, t))  &
                                      * moint2(ijkl2intindex(r, b, e, c))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(d, a, c, t))  &
                                      * moint2(ijkl2intindex(e, b, r, c))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(d, a, c, t))  &
                                      * moint2(ijkl2intindex(e, b, r, c))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(r)-lam(s)-lam(t)+lam(b)+lam(c)+lam(d)) *  &
                                    (-lam(s)-lam(t)+lam(d)+lam(e))  &
                                    )
                            end do
                        end do
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
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(r, a, s, b))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(r, a, s, b))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(s, a, r, b))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(s, a, r, b))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(d, t, b, u))  &
                                      * moint2(ijkl2intindex(r, a, s, b))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(d, t, b, u))  &
                                      * moint2(ijkl2intindex(s, a, r, b))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(c)) *  &
                                    (-lam(r)-lam(s)-lam(t)-lam(u)+lam(a)+lam(b)+lam(c)+lam(d)) *  &
                                    (-lam(t)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term4(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, a, d, t))  &
                                      * moint2(ijkl2intindex(r, b, s, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, a, d, t))  &
                                      * moint2(ijkl2intindex(r, b, s, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, a, d, t))  &
                                      * moint2(ijkl2intindex(s, b, r, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, a, d, t))  &
                                      * moint2(ijkl2intindex(s, b, r, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(d, a, b, t))  &
                                      * moint2(ijkl2intindex(r, b, s, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(d, a, b, t))  &
                                      * moint2(ijkl2intindex(r, b, s, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(d, a, b, t))  &
                                      * moint2(ijkl2intindex(s, b, r, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(d, a, b, t))  &
                                      * moint2(ijkl2intindex(s, b, r, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(c)) *  &
                                    (-lam(r)-lam(s)-lam(t)+lam(b)+lam(c)+lam(d)) *  &
                                    (-lam(t)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term5(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(r, b, s, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(r, b, s, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(s, b, r, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(s, b, r, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(r)-lam(s)-lam(t)+lam(b)+lam(c)+lam(d)) *  &
                                    (-lam(t)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term6(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, t, d, u))  &
                                      * moint2(ijkl2intindex(r, a, t, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, t, d, u))  &
                                      * moint2(ijkl2intindex(r, a, t, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, t, d, u))  &
                                      * moint2(ijkl2intindex(t, a, r, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, t, d, u))  &
                                      * moint2(ijkl2intindex(t, a, r, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(r)-lam(s)-lam(t)-lam(u)+lam(a)+lam(b)+lam(c)+lam(d)) *  &
                                    (-lam(s)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term7(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(t, a, d, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(t, a, d, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(d, a, t, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(d, a, t, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, t, r, u))  &
                                      * moint2(ijkl2intindex(t, a, d, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, t, r, u))  &
                                      * moint2(ijkl2intindex(t, a, d, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, t, r, u))  &
                                      * moint2(ijkl2intindex(d, a, t, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, t, r, u))  &
                                      * moint2(ijkl2intindex(d, a, t, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(s)-lam(t)-lam(u)+lam(a)+lam(b)+lam(c)) *  &
                                    (-lam(s)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term8(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(s, a, d, b))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(s, a, d, b))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(d, a, s, b))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(d, a, s, b))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, t, r, u))  &
                                      * moint2(ijkl2intindex(s, a, d, b))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, t, r, u))  &
                                      * moint2(ijkl2intindex(d, a, s, b))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(s)-lam(t)-lam(u)+lam(a)+lam(b)+lam(c)) *  &
                                    (-lam(t)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term9(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, r, s, t, u, v
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do r = Noccupied+1, size(lam)
                do s = Noccupied+1, size(lam)
                    do t = Noccupied+1, size(lam)
                        do u = Noccupied+1, size(lam)
                            do v = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(s, a, t, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(s, a, t, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(t, a, s, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(3+3) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(t, a, s, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, t, r, u))  &
                                      * moint2(ijkl2intindex(s, a, t, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, t, r, u))  &
                                      * moint2(ijkl2intindex(s, a, t, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, t, r, u))  &
                                      * moint2(ijkl2intindex(t, a, s, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, t, r, u))  &
                                      * moint2(ijkl2intindex(t, a, s, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(s, a, t, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(3+3) * 2**3  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(s, a, t, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(t, a, s, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(t, a, s, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(3+3) * 2**3  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, t, r, u))  &
                                      * moint2(ijkl2intindex(s, a, t, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, t, r, u))  &
                                      * moint2(ijkl2intindex(s, a, t, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, t, r, u))  &
                                      * moint2(ijkl2intindex(t, a, s, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, t, r, u))  &
                                      * moint2(ijkl2intindex(t, a, s, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(s)-lam(t)-lam(u)+lam(a)+lam(b)+lam(c)) *  &
                                    (-lam(u)-lam(v)+lam(b)+lam(c))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term10(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, e, r, s, t
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do e = 1, Noccupied
                    do r = Noccupied+1, size(lam)
                        do s = Noccupied+1, size(lam)
                            do t = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(t, b, e, c))  &
                                      * moint2(ijkl2intindex(r, d, s, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(t, b, e, c))  &
                                      * moint2(ijkl2intindex(s, d, r, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(e, b, t, c))  &
                                      * moint2(ijkl2intindex(r, d, s, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(e, b, t, c))  &
                                      * moint2(ijkl2intindex(s, d, r, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(d, a, c, t))  &
                                      * moint2(ijkl2intindex(t, b, e, c))  &
                                      * moint2(ijkl2intindex(r, d, s, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(d, a, c, t))  &
                                      * moint2(ijkl2intindex(t, b, e, c))  &
                                      * moint2(ijkl2intindex(s, d, r, e))  &
                                        + (-1)**(3+5) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(d, a, c, t))  &
                                      * moint2(ijkl2intindex(e, b, t, c))  &
                                      * moint2(ijkl2intindex(r, d, s, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(d, a, c, t))  &
                                      * moint2(ijkl2intindex(e, b, t, c))  &
                                      * moint2(ijkl2intindex(s, d, r, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(t, b, e, c))  &
                                      * moint2(ijkl2intindex(s, d, r, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(e, b, t, c))  &
                                      * moint2(ijkl2intindex(r, d, s, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(d, a, c, t))  &
                                      * moint2(ijkl2intindex(t, b, e, c))  &
                                      * moint2(ijkl2intindex(r, d, s, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(d, a, c, t))  &
                                      * moint2(ijkl2intindex(e, b, t, c))  &
                                      * moint2(ijkl2intindex(r, d, s, e))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(r)-lam(s)-lam(t)+lam(b)+lam(c)+lam(d)) *  &
                                    (-lam(r)-lam(s)+lam(d)+lam(e))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term11(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(4+4) * 2**4  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(r, a, t, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(r, a, t, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(t, a, r, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(t, a, r, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(d, t, b, u))  &
                                      * moint2(ijkl2intindex(r, a, t, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(d, t, b, u))  &
                                      * moint2(ijkl2intindex(r, a, t, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(d, t, b, u))  &
                                      * moint2(ijkl2intindex(t, a, r, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(d, t, b, u))  &
                                      * moint2(ijkl2intindex(t, a, r, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(r, a, t, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(r, a, t, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(t, a, r, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(t, a, r, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(d, t, b, u))  &
                                      * moint2(ijkl2intindex(r, a, t, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(d, t, b, u))  &
                                      * moint2(ijkl2intindex(r, a, t, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(d, t, b, u))  &
                                      * moint2(ijkl2intindex(t, a, r, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(d, t, b, u))  &
                                      * moint2(ijkl2intindex(t, a, r, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(c)) *  &
                                    (-lam(r)-lam(s)-lam(t)-lam(u)+lam(a)+lam(b)+lam(c)+lam(d)) *  &
                                    (-lam(s)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term12(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, a, d, t))  &
                                      * moint2(ijkl2intindex(r, b, t, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, a, d, t))  &
                                      * moint2(ijkl2intindex(r, b, t, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, a, d, t))  &
                                      * moint2(ijkl2intindex(t, b, r, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, a, d, t))  &
                                      * moint2(ijkl2intindex(t, b, r, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(d, a, b, t))  &
                                      * moint2(ijkl2intindex(r, b, t, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(d, a, b, t))  &
                                      * moint2(ijkl2intindex(r, b, t, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(d, a, b, t))  &
                                      * moint2(ijkl2intindex(t, b, r, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(d, a, b, t))  &
                                      * moint2(ijkl2intindex(t, b, r, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, d, t))  &
                                      * moint2(ijkl2intindex(r, b, t, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, d, t))  &
                                      * moint2(ijkl2intindex(r, b, t, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, d, t))  &
                                      * moint2(ijkl2intindex(t, b, r, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, d, t))  &
                                      * moint2(ijkl2intindex(t, b, r, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(d, a, b, t))  &
                                      * moint2(ijkl2intindex(r, b, t, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(d, a, b, t))  &
                                      * moint2(ijkl2intindex(r, b, t, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(d, a, b, t))  &
                                      * moint2(ijkl2intindex(t, b, r, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(d, a, b, t))  &
                                      * moint2(ijkl2intindex(t, b, r, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(c)) *  &
                                    (-lam(r)-lam(s)-lam(t)+lam(b)+lam(c)+lam(d)) *  &
                                    (-lam(s)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term13(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(r, b, t, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(r, b, t, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(t, b, r, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(t, b, r, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(r, b, t, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(r, b, t, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(t, b, r, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, d, t))  &
                                      * moint2(ijkl2intindex(t, b, r, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(r)-lam(s)-lam(t)+lam(b)+lam(c)+lam(d)) *  &
                                    (-lam(s)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term14(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(r, t, b, u))  &
                                      * moint2(ijkl2intindex(t, a, d, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(r, t, b, u))  &
                                      * moint2(ijkl2intindex(t, a, d, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(r, t, b, u))  &
                                      * moint2(ijkl2intindex(d, a, t, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(r, t, b, u))  &
                                      * moint2(ijkl2intindex(d, a, t, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, r, u))  &
                                      * moint2(ijkl2intindex(t, a, d, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, r, u))  &
                                      * moint2(ijkl2intindex(t, a, d, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, r, u))  &
                                      * moint2(ijkl2intindex(d, a, t, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, r, u))  &
                                      * moint2(ijkl2intindex(d, a, t, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(r, t, b, u))  &
                                      * moint2(ijkl2intindex(t, a, d, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(r, t, b, u))  &
                                      * moint2(ijkl2intindex(t, a, d, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(r, t, b, u))  &
                                      * moint2(ijkl2intindex(d, a, t, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(r, t, b, u))  &
                                      * moint2(ijkl2intindex(d, a, t, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, t, r, u))  &
                                      * moint2(ijkl2intindex(t, a, d, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, t, r, u))  &
                                      * moint2(ijkl2intindex(t, a, d, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, t, r, u))  &
                                      * moint2(ijkl2intindex(d, a, t, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, t, r, u))  &
                                      * moint2(ijkl2intindex(d, a, t, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(c)) *  &
                                    (-lam(s)-lam(t)-lam(u)+lam(a)+lam(b)+lam(c)) *  &
                                    (-lam(s)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term15(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, e, r, s, t
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do e = 1, Noccupied
                    do r = Noccupied+1, size(lam)
                        do s = Noccupied+1, size(lam)
                            do t = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(3+5) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(d, b, e, c))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(d, b, e, c))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, e, c))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, e, c))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(d, b, e, c))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(d, b, e, c))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, e, c))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, e, c))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(s)-lam(t)+lam(b)+lam(c)) *  &
                                    (-lam(s)-lam(t)+lam(d)+lam(e))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term16(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(4+4) * 2**4  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(r, a, b, t))  &
                                      * moint2(ijkl2intindex(t, b, d, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(r, a, b, t))  &
                                      * moint2(ijkl2intindex(t, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(r, a, b, t))  &
                                      * moint2(ijkl2intindex(d, b, t, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(r, a, b, t))  &
                                      * moint2(ijkl2intindex(d, b, t, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, a, r, t))  &
                                      * moint2(ijkl2intindex(t, b, d, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, a, r, t))  &
                                      * moint2(ijkl2intindex(t, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, t, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, t, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, b, t))  &
                                      * moint2(ijkl2intindex(t, b, d, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, b, t))  &
                                      * moint2(ijkl2intindex(t, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, b, t))  &
                                      * moint2(ijkl2intindex(d, b, t, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, b, t))  &
                                      * moint2(ijkl2intindex(d, b, t, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, r, t))  &
                                      * moint2(ijkl2intindex(t, b, d, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, r, t))  &
                                      * moint2(ijkl2intindex(t, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, t, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, t, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(c)) *  &
                                    (-lam(s)-lam(t)+lam(b)+lam(c)) *  &
                                    (-lam(s)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term17(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(t, b, d, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(t, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(d, b, t, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(d, b, t, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(t, b, d, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(t, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, t, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, t, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(t, b, d, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(t, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(d, b, t, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(d, b, t, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(t, b, d, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(t, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, t, u))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, t, u))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(s)-lam(t)+lam(b)+lam(c)) *  &
                                    (-lam(s)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term18(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(r, t, b, u))  &
                                      * moint2(ijkl2intindex(s, a, d, b))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(r, t, b, u))  &
                                      * moint2(ijkl2intindex(s, a, d, b))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(r, t, b, u))  &
                                      * moint2(ijkl2intindex(d, a, s, b))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(r, t, b, u))  &
                                      * moint2(ijkl2intindex(d, a, s, b))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, r, u))  &
                                      * moint2(ijkl2intindex(s, a, d, b))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, r, u))  &
                                      * moint2(ijkl2intindex(d, a, s, b))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(r, t, b, u))  &
                                      * moint2(ijkl2intindex(s, a, d, b))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(r, t, b, u))  &
                                      * moint2(ijkl2intindex(s, a, d, b))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(r, t, b, u))  &
                                      * moint2(ijkl2intindex(d, a, s, b))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(r, t, b, u))  &
                                      * moint2(ijkl2intindex(d, a, s, b))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, t, r, u))  &
                                      * moint2(ijkl2intindex(s, a, d, b))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, t, r, u))  &
                                      * moint2(ijkl2intindex(d, a, s, b))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(c)) *  &
                                    (-lam(s)-lam(t)-lam(u)+lam(a)+lam(b)+lam(c)) *  &
                                    (-lam(t)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term19(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, r, s, t, u, v
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do r = Noccupied+1, size(lam)
                do s = Noccupied+1, size(lam)
                    do t = Noccupied+1, size(lam)
                        do u = Noccupied+1, size(lam)
                            do v = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, c, s))  &
                                      * moint2(ijkl2intindex(r, t, a, u))  &
                                      * moint2(ijkl2intindex(s, a, t, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, c, s))  &
                                      * moint2(ijkl2intindex(r, t, a, u))  &
                                      * moint2(ijkl2intindex(s, a, t, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, c, s))  &
                                      * moint2(ijkl2intindex(r, t, a, u))  &
                                      * moint2(ijkl2intindex(t, a, s, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, c, s))  &
                                      * moint2(ijkl2intindex(r, t, a, u))  &
                                      * moint2(ijkl2intindex(t, a, s, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, c, s))  &
                                      * moint2(ijkl2intindex(a, t, r, u))  &
                                      * moint2(ijkl2intindex(s, a, t, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, c, s))  &
                                      * moint2(ijkl2intindex(a, t, r, u))  &
                                      * moint2(ijkl2intindex(s, a, t, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(3+3) * 2**3  &
                                      * moint2(ijkl2intindex(b, r, c, s))  &
                                      * moint2(ijkl2intindex(a, t, r, u))  &
                                      * moint2(ijkl2intindex(t, a, s, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, c, s))  &
                                      * moint2(ijkl2intindex(a, t, r, u))  &
                                      * moint2(ijkl2intindex(t, a, s, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(b)+lam(c)) *  &
                                    (-lam(s)-lam(t)-lam(u)+lam(a)+lam(b)+lam(c)) *  &
                                    (-lam(u)-lam(v)+lam(b)+lam(c))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term20(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(r, a, b, t))  &
                                      * moint2(ijkl2intindex(s, b, d, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(r, a, b, t))  &
                                      * moint2(ijkl2intindex(s, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(r, a, b, t))  &
                                      * moint2(ijkl2intindex(d, b, s, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(r, a, b, t))  &
                                      * moint2(ijkl2intindex(d, b, s, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, a, r, t))  &
                                      * moint2(ijkl2intindex(s, b, d, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, a, r, t))  &
                                      * moint2(ijkl2intindex(s, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, s, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, s, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, b, t))  &
                                      * moint2(ijkl2intindex(s, b, d, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, b, t))  &
                                      * moint2(ijkl2intindex(s, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, b, t))  &
                                      * moint2(ijkl2intindex(d, b, s, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, b, t))  &
                                      * moint2(ijkl2intindex(d, b, s, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, r, t))  &
                                      * moint2(ijkl2intindex(s, b, d, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, r, t))  &
                                      * moint2(ijkl2intindex(s, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, s, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, r, t))  &
                                      * moint2(ijkl2intindex(d, b, s, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(c)) *  &
                                    (-lam(s)-lam(t)+lam(b)+lam(c)) *  &
                                    (-lam(t)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term21(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, r, s, t, u, v
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do r = Noccupied+1, size(lam)
                do s = Noccupied+1, size(lam)
                    do t = Noccupied+1, size(lam)
                        do u = Noccupied+1, size(lam)
                            do v = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(3+3) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(s, u, t, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(s, u, t, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(t, u, s, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(s, u, t, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(s, u, t, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(t, u, s, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(s, u, t, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(s, u, t, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, a, c, t))  &
                                      * moint2(ijkl2intindex(t, u, s, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(s, u, t, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(s, u, t, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(c, a, r, t))  &
                                      * moint2(ijkl2intindex(t, u, s, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(s)-lam(t)+lam(b)+lam(c)) *  &
                                    (-lam(u)-lam(v)+lam(b)+lam(c))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term22(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, e, r, s, t
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do e = 1, Noccupied
                    do r = Noccupied+1, size(lam)
                        do s = Noccupied+1, size(lam)
                            do t = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, d, s))  &
                                      * moint2(ijkl2intindex(b, a, c, t))  &
                                      * moint2(ijkl2intindex(t, b, e, c))  &
                                      * moint2(ijkl2intindex(r, d, s, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, d, s))  &
                                      * moint2(ijkl2intindex(b, a, c, t))  &
                                      * moint2(ijkl2intindex(t, b, e, c))  &
                                      * moint2(ijkl2intindex(s, d, r, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, d, s))  &
                                      * moint2(ijkl2intindex(b, a, c, t))  &
                                      * moint2(ijkl2intindex(e, b, t, c))  &
                                      * moint2(ijkl2intindex(r, d, s, e))  &
                                        + (-1)**(3+5) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, d, s))  &
                                      * moint2(ijkl2intindex(b, a, c, t))  &
                                      * moint2(ijkl2intindex(e, b, t, c))  &
                                      * moint2(ijkl2intindex(s, d, r, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(d, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, c, t))  &
                                      * moint2(ijkl2intindex(t, b, e, c))  &
                                      * moint2(ijkl2intindex(s, d, r, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(d, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, c, t))  &
                                      * moint2(ijkl2intindex(e, b, t, c))  &
                                      * moint2(ijkl2intindex(s, d, r, e))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(d)) *  &
                                    (-lam(r)-lam(s)-lam(t)+lam(b)+lam(c)+lam(d)) *  &
                                    (-lam(r)-lam(s)+lam(d)+lam(e))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term23(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, d, s))  &
                                      * moint2(ijkl2intindex(a, t, b, u))  &
                                      * moint2(ijkl2intindex(r, a, t, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, d, s))  &
                                      * moint2(ijkl2intindex(a, t, b, u))  &
                                      * moint2(ijkl2intindex(r, a, t, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, d, s))  &
                                      * moint2(ijkl2intindex(a, t, b, u))  &
                                      * moint2(ijkl2intindex(t, a, r, b))  &
                                      * moint2(ijkl2intindex(s, c, u, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(c, r, d, s))  &
                                      * moint2(ijkl2intindex(a, t, b, u))  &
                                      * moint2(ijkl2intindex(t, a, r, b))  &
                                      * moint2(ijkl2intindex(u, c, s, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(c)+lam(d)) *  &
                                    (-lam(r)-lam(s)-lam(t)-lam(u)+lam(a)+lam(b)+lam(c)+lam(d)) *  &
                                    (-lam(s)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term24(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, e, r, s, t
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do e = 1, Noccupied
                    do r = Noccupied+1, size(lam)
                        do s = Noccupied+1, size(lam)
                            do t = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(3+5) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, d, s))  &
                                      * moint2(ijkl2intindex(b, a, c, t))  &
                                      * moint2(ijkl2intindex(r, b, e, c))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, d, s))  &
                                      * moint2(ijkl2intindex(b, a, c, t))  &
                                      * moint2(ijkl2intindex(r, b, e, c))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, d, s))  &
                                      * moint2(ijkl2intindex(b, a, c, t))  &
                                      * moint2(ijkl2intindex(e, b, r, c))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, d, s))  &
                                      * moint2(ijkl2intindex(b, a, c, t))  &
                                      * moint2(ijkl2intindex(e, b, r, c))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(d, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, c, t))  &
                                      * moint2(ijkl2intindex(r, b, e, c))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(d, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, c, t))  &
                                      * moint2(ijkl2intindex(r, b, e, c))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(d, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, c, t))  &
                                      * moint2(ijkl2intindex(e, b, r, c))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(d, r, a, s))  &
                                      * moint2(ijkl2intindex(b, a, c, t))  &
                                      * moint2(ijkl2intindex(e, b, r, c))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(d)) *  &
                                    (-lam(r)-lam(s)-lam(t)+lam(b)+lam(c)+lam(d)) *  &
                                    (-lam(s)-lam(t)+lam(d)+lam(e))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term25(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(t, a, u, b))  &
                                      * moint2(ijkl2intindex(r, c, s, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(t, a, u, b))  &
                                      * moint2(ijkl2intindex(s, c, r, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(u, a, t, b))  &
                                      * moint2(ijkl2intindex(r, c, s, d))  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, c, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(u, a, t, b))  &
                                      * moint2(ijkl2intindex(s, c, r, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(t, a, u, b))  &
                                      * moint2(ijkl2intindex(s, c, r, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(c, r, a, s))  &
                                      * moint2(ijkl2intindex(b, t, d, u))  &
                                      * moint2(ijkl2intindex(u, a, t, b))  &
                                      * moint2(ijkl2intindex(s, c, r, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(c)) *  &
                                    (-lam(r)-lam(s)-lam(t)-lam(u)+lam(a)+lam(b)+lam(c)+lam(d)) *  &
                                    (-lam(r)-lam(s)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term26(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, r, s, t, u, v
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do r = Noccupied+1, size(lam)
                do s = Noccupied+1, size(lam)
                    do t = Noccupied+1, size(lam)
                        do u = Noccupied+1, size(lam)
                            do v = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, c, s))  &
                                      * moint2(ijkl2intindex(r, t, a, u))  &
                                      * moint2(ijkl2intindex(t, a, u, v))  &
                                      * moint2(ijkl2intindex(s, b, v, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, c, s))  &
                                      * moint2(ijkl2intindex(r, t, a, u))  &
                                      * moint2(ijkl2intindex(t, a, u, v))  &
                                      * moint2(ijkl2intindex(v, b, s, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, c, s))  &
                                      * moint2(ijkl2intindex(r, t, a, u))  &
                                      * moint2(ijkl2intindex(u, a, t, v))  &
                                      * moint2(ijkl2intindex(s, b, v, c))  &
                                        + (-1)**(3+3) * 2**3  &
                                      * moint2(ijkl2intindex(b, r, c, s))  &
                                      * moint2(ijkl2intindex(r, t, a, u))  &
                                      * moint2(ijkl2intindex(u, a, t, v))  &
                                      * moint2(ijkl2intindex(v, b, s, c))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(b)+lam(c)) *  &
                                    (-lam(s)-lam(t)-lam(u)+lam(a)+lam(b)+lam(c)) *  &
                                    (-lam(s)-lam(v)+lam(b)+lam(c))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term27(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, r, s, t, u, v
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do r = Noccupied+1, size(lam)
                do s = Noccupied+1, size(lam)
                    do t = Noccupied+1, size(lam)
                        do u = Noccupied+1, size(lam)
                            do v = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(3+3) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(t, a, u, v))  &
                                      * moint2(ijkl2intindex(s, b, v, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(t, a, u, v))  &
                                      * moint2(ijkl2intindex(v, b, s, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(u, a, t, v))  &
                                      * moint2(ijkl2intindex(s, b, v, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(u, a, t, v))  &
                                      * moint2(ijkl2intindex(v, b, s, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(t, a, u, v))  &
                                      * moint2(ijkl2intindex(s, b, v, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(t, a, u, v))  &
                                      * moint2(ijkl2intindex(v, b, s, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(u, a, t, v))  &
                                      * moint2(ijkl2intindex(s, b, v, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(b, r, a, s))  &
                                      * moint2(ijkl2intindex(r, t, c, u))  &
                                      * moint2(ijkl2intindex(u, a, t, v))  &
                                      * moint2(ijkl2intindex(v, b, s, c))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(s)-lam(t)-lam(u)+lam(a)+lam(b)+lam(c)) *  &
                                    (-lam(s)-lam(v)+lam(b)+lam(c))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term28(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, e, r, s, t
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do e = 1, Noccupied
                    do r = Noccupied+1, size(lam)
                        do s = Noccupied+1, size(lam)
                            do t = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(3+5) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, b))  &
                                      * moint2(ijkl2intindex(r, c, e, t))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, b))  &
                                      * moint2(ijkl2intindex(r, c, e, t))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, b))  &
                                      * moint2(ijkl2intindex(e, c, r, t))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, b))  &
                                      * moint2(ijkl2intindex(e, c, r, t))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(d, a, c, b))  &
                                      * moint2(ijkl2intindex(r, c, e, t))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(d, a, c, b))  &
                                      * moint2(ijkl2intindex(r, c, e, t))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(d, a, c, b))  &
                                      * moint2(ijkl2intindex(e, c, r, t))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(d, a, c, b))  &
                                      * moint2(ijkl2intindex(e, c, r, t))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(r)-lam(s)+lam(c)+lam(d)) *  &
                                    (-lam(s)-lam(t)+lam(d)+lam(e))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term29(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, e, r, s, t
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do e = 1, Noccupied
                    do r = Noccupied+1, size(lam)
                        do s = Noccupied+1, size(lam)
                            do t = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(3+5) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, b))  &
                                      * moint2(ijkl2intindex(d, c, e, t))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, b))  &
                                      * moint2(ijkl2intindex(d, c, e, t))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                        + (-1)**(2+5) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, b))  &
                                      * moint2(ijkl2intindex(d, c, e, t))  &
                                      * moint2(ijkl2intindex(s, d, t, e))  &
                                        + (-1)**(1+5) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, b))  &
                                      * moint2(ijkl2intindex(d, c, e, t))  &
                                      * moint2(ijkl2intindex(t, d, s, e))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(s)+lam(c)) *  &
                                    (-lam(s)-lam(t)+lam(d)+lam(e))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term30(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, b))  &
                                      * moint2(ijkl2intindex(s, t, d, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, b))  &
                                      * moint2(ijkl2intindex(s, t, d, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, c, b))  &
                                      * moint2(ijkl2intindex(d, t, s, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, b))  &
                                      * moint2(ijkl2intindex(s, t, d, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, b))  &
                                      * moint2(ijkl2intindex(s, t, d, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, r, b))  &
                                      * moint2(ijkl2intindex(d, t, s, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(s)+lam(c)) *  &
                                    (-lam(t)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term31(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, e, f, r, s
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do e = 1, Noccupied
                    do f = 1, Noccupied
                        do r = Noccupied+1, size(lam)
                            do s = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(2+6) * 2**2 / 2._dp  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, b))  &
                                      * moint2(ijkl2intindex(e, c, f, d))  &
                                      * moint2(ijkl2intindex(r, e, s, f))  &
                                        + (-1)**(1+6) * 2**1 / 2._dp  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, b))  &
                                      * moint2(ijkl2intindex(e, c, f, d))  &
                                      * moint2(ijkl2intindex(s, e, r, f))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(r)-lam(s)+lam(c)+lam(d)) *  &
                                    (-lam(r)-lam(s)+lam(e)+lam(f))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term32(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(2+4) * 2**2 / 2._dp  &
                                      * moint2(ijkl2intindex(c, r, d, s))  &
                                      * moint2(ijkl2intindex(a, t, b, u))  &
                                      * moint2(ijkl2intindex(r, a, s, b))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1 / 2._dp  &
                                      * moint2(ijkl2intindex(c, r, d, s))  &
                                      * moint2(ijkl2intindex(a, t, b, u))  &
                                      * moint2(ijkl2intindex(r, a, s, b))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(c)+lam(d)) *  &
                                    (-lam(r)-lam(s)-lam(t)-lam(u)+lam(a)+lam(b)+lam(c)+lam(d)) *  &
                                    (-lam(t)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term33(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(2+4) * 2**2 / 2._dp  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, b))  &
                                      * moint2(ijkl2intindex(r, t, s, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1 / 2._dp  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, a, d, b))  &
                                      * moint2(ijkl2intindex(r, t, s, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(r)-lam(s)+lam(c)+lam(d)) *  &
                                    (-lam(t)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term34(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, r, s, t, u, v
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do r = Noccupied+1, size(lam)
                do s = Noccupied+1, size(lam)
                    do t = Noccupied+1, size(lam)
                        do u = Noccupied+1, size(lam)
                            do v = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(3+3) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, s, u))  &
                                      * moint2(ijkl2intindex(t, a, c, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, s, u))  &
                                      * moint2(ijkl2intindex(t, a, c, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, s, u))  &
                                      * moint2(ijkl2intindex(c, a, t, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, s, u))  &
                                      * moint2(ijkl2intindex(c, a, t, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(s, t, r, u))  &
                                      * moint2(ijkl2intindex(t, a, c, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(s, t, r, u))  &
                                      * moint2(ijkl2intindex(t, a, c, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(s, t, r, u))  &
                                      * moint2(ijkl2intindex(c, a, t, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(s, t, r, u))  &
                                      * moint2(ijkl2intindex(c, a, t, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(t)-lam(u)+lam(a)+lam(b)) *  &
                                    (-lam(u)-lam(v)+lam(b)+lam(c))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term35(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(3+4) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, s, t))  &
                                      * moint2(ijkl2intindex(c, b, d, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, s, t))  &
                                      * moint2(ijkl2intindex(c, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                        + (-1)**(2+4) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(s, a, r, t))  &
                                      * moint2(ijkl2intindex(c, b, d, u))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(s, a, r, t))  &
                                      * moint2(ijkl2intindex(c, b, d, u))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(t)+lam(b)) *  &
                                    (-lam(t)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term36(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, r, s, t, u, v
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do r = Noccupied+1, size(lam)
                do s = Noccupied+1, size(lam)
                    do t = Noccupied+1, size(lam)
                        do u = Noccupied+1, size(lam)
                            do v = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(3+3) * 2**3  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, s, t))  &
                                      * moint2(ijkl2intindex(t, u, c, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, s, t))  &
                                      * moint2(ijkl2intindex(t, u, c, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, a, s, t))  &
                                      * moint2(ijkl2intindex(c, u, t, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(2+3) * 2**2  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(s, a, r, t))  &
                                      * moint2(ijkl2intindex(t, u, c, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(s, a, r, t))  &
                                      * moint2(ijkl2intindex(t, u, c, v))  &
                                      * moint2(ijkl2intindex(v, b, u, c))  &
                                        + (-1)**(1+3) * 2**1  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(s, a, r, t))  &
                                      * moint2(ijkl2intindex(c, u, t, v))  &
                                      * moint2(ijkl2intindex(u, b, v, c))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(t)+lam(b)) *  &
                                    (-lam(u)-lam(v)+lam(b)+lam(c))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term37(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(2+4) * 2**2 / 2._dp  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, t, d, u))  &
                                      * moint2(ijkl2intindex(t, a, u, b))  &
                                      * moint2(ijkl2intindex(r, c, s, d))  &
                                        + (-1)**(1+4) * 2**1 / 2._dp  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(c, t, d, u))  &
                                      * moint2(ijkl2intindex(t, a, u, b))  &
                                      * moint2(ijkl2intindex(s, c, r, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(r)-lam(s)-lam(t)-lam(u)+lam(a)+lam(b)+lam(c)+lam(d)) *  &
                                    (-lam(r)-lam(s)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term38(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, c, d, r, s, t, u
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do c = 1, Noccupied
            do d = 1, Noccupied
                do r = Noccupied+1, size(lam)
                    do s = Noccupied+1, size(lam)
                        do t = Noccupied+1, size(lam)
                            do u = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(2+4) * 2**2 / 2._dp  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, s, u))  &
                                      * moint2(ijkl2intindex(c, a, d, b))  &
                                      * moint2(ijkl2intindex(t, c, u, d))  &
                                        + (-1)**(1+4) * 2**1 / 2._dp  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, s, u))  &
                                      * moint2(ijkl2intindex(c, a, d, b))  &
                                      * moint2(ijkl2intindex(u, c, t, d))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(t)-lam(u)+lam(a)+lam(b)) *  &
                                    (-lam(t)-lam(u)+lam(c)+lam(d))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

real(dp) function term39(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
integer :: a, b, r, s, t, u, v, w
E2 = 0
do a = 1, Noccupied
    do b = 1, Noccupied
        do r = Noccupied+1, size(lam)
            do s = Noccupied+1, size(lam)
                do t = Noccupied+1, size(lam)
                    do u = Noccupied+1, size(lam)
                        do v = Noccupied+1, size(lam)
                            do w = Noccupied+1, size(lam)
                                E2 = E2 + (  &
                                        + (-1)**(2+2) * 2**2 / 2._dp  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, s, u))  &
                                      * moint2(ijkl2intindex(t, v, u, w))  &
                                      * moint2(ijkl2intindex(v, a, w, b))  &
                                        + (-1)**(1+2) * 2**1 / 2._dp  &
                                      * moint2(ijkl2intindex(a, r, b, s))  &
                                      * moint2(ijkl2intindex(r, t, s, u))  &
                                      * moint2(ijkl2intindex(t, v, u, w))  &
                                      * moint2(ijkl2intindex(w, a, v, b))  &
                                    ) / (  &
                                    (-lam(r)-lam(s)+lam(a)+lam(b)) *  &
                                    (-lam(t)-lam(u)+lam(a)+lam(b)) *  &
                                    (-lam(v)-lam(w)+lam(a)+lam(b))  &
                                    )
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end do
end function

end module
