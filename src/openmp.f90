module openmp
use types, only: dp
implicit none
private
public omp_get_thread_num, omp_get_num_threads, omp_get_max_threads, &
    omp_get_wtime

interface

    integer function omp_get_thread_num()
    end function

    integer function omp_get_num_threads()
    end function

    integer function omp_get_max_threads()
    end function

    real(dp) function omp_get_wtime()
    import :: dp
    end function

end interface

end module
