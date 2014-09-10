module openmp
use types, only: dp
implicit none
private
public omp_get_thread_num, omp_get_wtime

contains

    integer function omp_get_thread_num()
    omp_get_thread_num = 0
    end function

    real(dp) function omp_get_wtime()
    omp_get_wtime = 0
    end function

end module
