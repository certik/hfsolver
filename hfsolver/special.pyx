cimport c_hfsolver

def Inu(int k, double x):
    cdef double r
    c_hfsolver.hfsolver_Inu(&k, &x, &r)
    return r

def Knu(int k, double x):
    cdef double r
    c_hfsolver.hfsolver_Knu(&k, &x, &r)
    return r
