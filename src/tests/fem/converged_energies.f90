module converged_energies

! This module stores energies for various test problems that both the FE and FFT
! solvers must converged to (the solvers are in different programs that share
! this module).

use types, only: dp
implicit none

! [Ts, Een, Eee, Exc]
! nuclear: alpha = 6, electronic alpha = 5
real(dp), parameter :: one_gaussian(*) = [10.61904507_dp, -2.92172113_dp, &
            1.30109500_dp, -1.43805890_dp]

! nuclear: alpha = 6, electronic alpha = 5
real(dp), parameter :: four_gaussians(*) = [107.03329544_dp, 1.17805865_dp, &
            20.81751993_dp, -8.98245221_dp]

end module
