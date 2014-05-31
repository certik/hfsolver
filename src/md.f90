module md
use types, only: dp
implicit none
private
public velocity_verlet, minimize_energy, unfold_positions

interface
    subroutine forces_func(X, f)
    ! Calculate forces for particles at positions X
    import :: dp
    implicit none
    ! (i, j) is the i-th particle, j-component (j=1, 2, 3)
    real(dp), intent(in) :: X(:, :) ! positions
    real(dp), intent(out) :: f(:, :) ! forces
    end subroutine

    real(dp) function energy_func(X)
    ! Calculate the energy for particles at positions X
    import :: dp
    implicit none
    ! (i, j) is the i-th particle, j-component (j=1, 2, 3)
    real(dp), intent(in) :: X(:, :) ! positions
    end function
end interface

contains

subroutine velocity_verlet(dt, m, forces, f, V, X)
! Propagates f, V, X: t -> t+dt
! On input, f, V, X are given at time t, on output they are given at time t+dt.
! X(i, j) is i-th particle, j comp. (j=1, 2, 3)
real(dp), intent(in) :: dt ! time step
real(dp), intent(in) :: m(:) ! m(i) the mass of i-th particle
! Callback that calculates forces from positions
procedure(forces_func) :: forces
real(dp), intent(inout) :: f(:, :) ! f(t) -> f(t+dt)
real(dp), intent(inout) :: V(:, :) ! V(t) -> V(t+dt)
real(dp), intent(inout) :: X(:, :) ! X(t) -> X(t+dt)
real(dp) :: f_next(3, size(f, 2))
X = X + V*dt + f*dt**2/(2*spread(m, 1, 3))
call forces(X, f_next)
V = V + (f + f_next)*dt/(2*spread(m, 1, 3))
f = f_next
end subroutine

subroutine minimize_energy(forces, energy, X, f, h0, max_iter)
! Callback that calculates forces from positions
procedure(forces_func) :: forces
procedure(energy_func) :: energy
real(dp), intent(inout) :: X(:, :)
real(dp), intent(out) :: f(:, :) ! The final force 'f'
real(dp), intent(in) :: h0
integer, intent(in) :: max_iter
real(dp) :: h, E, Enew, Xnew(3, size(X, 2))
integer :: i
call forces(X, f)
h = h0
E = energy(X)
call forces(X, f)
do i = 1, max_iter
    print *, i, E, h
    Xnew = X + f/maxval(sqrt(sum(f**2, dim=1))) * h
    Enew = energy(Xnew)
    if (Enew < E) then
        ! new positions are accepted
        X = Xnew
        E = Enew
        call forces(X, f)
        h = 1.2_dp * h
    else
        ! new positions are rejected
        h = 0.2_dp * h
    end if
end do
end subroutine

subroutine unfold_positions(L, X, Xu)
! Unwinds periodic positions. It is using the positions of particles from the
! first time step as a reference and then tracks them as they evolve possibly
! outside of this box. The X and Xu arrays are of the type X(:, j, i), which
! are the (x,y,z) coordinates of the j-th particle in the i-th time step.
real(dp), intent(in) :: L ! Box length
! Positions in [0, L]^3 with possible jumps:
real(dp), intent(in) :: X(:, :, :)
! Unwinded positions in (-oo, oo)^3 with no jumps (continuous):
real(dp), intent(out) :: Xu(:, :, :)
real(dp) :: d(3), Xj(3)
integer :: i, j
Xu(:, :, 1) = X(:, :, 1)
do i = 2, size(X, 3)
    do j = 1, size(X, 2)
        Xj = X(:, j, i-1) - X(:, j, i) + [L/2, L/2, L/2]
        Xj = Xj - L*floor(Xj/L)
        d = [L/2, L/2, L/2] - Xj
        Xu(:, j, i) = Xu(:, j, i-1) + d
    end do
end do
end subroutine

end module
