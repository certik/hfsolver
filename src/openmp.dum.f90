module openmp
implicit none
private
public omp_get_thread_num

contains

    integer function omp_get_thread_num()
    omp_get_thread_num = 0
    end function

end module
