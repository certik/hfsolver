module openmp
implicit none
private
public omp_get_thread_num

interface

    integer function omp_get_thread_num()
    end function

end interface

end module
