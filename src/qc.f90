module qc
implicit none


interface

    subroutine getS(nbf, nprim, istart, xcenter, ycenter, zcenter, &
        lpower, mpower, npower, n2, coef, alpha, S) bind(c, name="getS")
    use iso_c_binding, only: c_int, c_double
    implicit none
    integer(c_int), value, intent(in) :: nbf, n2
    integer(c_int), intent(in) :: nprim(nbf)
    integer(c_int), intent(in) :: istart(nbf)
    real(c_double), intent(in) :: xcenter(nbf), ycenter(nbf), zcenter(nbf)
    integer(c_int), intent(in) :: lpower(nbf), mpower(nbf), npower(nbf)
    real(c_double), dimension(n2), intent(in) :: coef, alpha
    real(c_double), intent(out) :: S(nbf**2)
    end subroutine

    subroutine getT(nbf, nprim, istart, xcenter, ycenter, zcenter, &
        lpower, mpower, npower, n2, coef, alpha, T) bind(c, name="getT")
    use iso_c_binding, only: c_int, c_double
    implicit none
    integer(c_int), value, intent(in) :: nbf, n2
    integer(c_int), intent(in) :: nprim(nbf)
    integer(c_int), intent(in) :: istart(nbf)
    real(c_double), intent(in) :: xcenter(nbf), ycenter(nbf), zcenter(nbf)
    integer(c_int), intent(in) :: lpower(nbf), mpower(nbf), npower(nbf)
    real(c_double), dimension(n2), intent(in) :: coef, alpha
    real(c_double), intent(out) :: T(nbf**2)
    end subroutine

    subroutine getV(nbf, nprim, istart, xcenter, ycenter, zcenter, &
        lpower, mpower, npower, n2, coef, alpha, &
        nat, atno, x, y, z, &
        V) bind(c, name="getV")
    use iso_c_binding, only: c_int, c_double
    implicit none
    integer(c_int), value, intent(in) :: nbf, n2
    integer(c_int), intent(in) :: nprim(nbf)
    integer(c_int), intent(in) :: istart(nbf)
    real(c_double), intent(in) :: xcenter(nbf), ycenter(nbf), zcenter(nbf)
    integer(c_int), intent(in) :: lpower(nbf), mpower(nbf), npower(nbf)
    real(c_double), dimension(n2), intent(in) :: coef, alpha
    integer(c_int), value, intent(in) :: nat
    integer(c_int), intent(in) :: atno(nat)
    real(c_double), intent(in) :: x(nat), y(nat), z(nat)
    real(c_double), intent(out) :: V(nbf**2)
    end subroutine

    subroutine getInts(nbf, nprim, istart, xcenter, ycenter, zcenter, &
        lpower, mpower, npower, n2, coef, alpha, &
        Ints) bind(c,name="getInts")
    use iso_c_binding, only: c_int, c_double
    implicit none
    integer(c_int), value, intent(in) :: nbf, n2
    integer(c_int), intent(in) :: nprim(nbf)
    integer(c_int), intent(in) :: istart(nbf)
    real(c_double), intent(in) :: xcenter(nbf), ycenter(nbf), zcenter(nbf)
    integer(c_int), intent(in) :: lpower(nbf), mpower(nbf), npower(nbf)
    real(c_double), dimension(n2), intent(in) :: coef, alpha
    ! The number of integrals in Ints is m*(m+1)/2 where m = nbf*(nbf+1)/2
    real(c_double), intent(out) :: Ints(nbf*(nbf**3 + 2*nbf**2 + 3*nbf + 2)/8)
    end subroutine

    function coulomb_repulsion( &
            xa, ya, za, la, ma, na, alphaa, &
            xb, yb, zb, lb, mb, nb, alphab, &
            xc, yc, zc, lc, mc, nc, alphac, &
            xd, yd, zd, ld, md, nd, alphad) result(r) bind(c)
    use iso_c_binding, only: c_int, c_double
    implicit none
    real(c_double), value, intent(in) :: xa, ya, za, alphaa
    integer(c_int), value, intent(in) :: la, ma, na
    real(c_double), value, intent(in) :: xb, yb, zb, alphab
    integer(c_int), value, intent(in) :: lb, mb, nb
    real(c_double), value, intent(in) :: xc, yc, zc, alphac
    integer(c_int), value, intent(in) :: lc, mc, nc
    real(c_double), value, intent(in) :: xd, yd, zd, alphad
    integer(c_int), value, intent(in) :: ld, md, nd
    real(c_double) :: r
    end function

end interface

end module
