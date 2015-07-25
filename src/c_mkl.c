#include <complex.h>

#define MKL_Complex16 double complex

#include "mkl_dfti.h"

/*
 * sizeof(x) == n1*n2*n3*sizeof(double complex)
 * This assumes Fortran style column wise ordering.
 */
int fft3_inplace(double complex *x, int n1, int n2, int n3)
{
    MKL_LONG status = 0;
    DFTI_DESCRIPTOR_HANDLE hand = 0;
    {
        MKL_LONG N[3]; N[0] = n3; N[1] = n2; N[2] = n1;
        status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_COMPLEX, 3, N);
        if (0 != status) {
            DftiFreeDescriptor(&hand);
            return 1;
        }
    }
    status = DftiCommitDescriptor(hand);
    if (0 != status) {
        DftiFreeDescriptor(&hand);
        return 1;
    }
    MKL_Complex16 *x_ = x;
    status = DftiComputeForward(hand, x_);
    if (0 != status) {
        DftiFreeDescriptor(&hand);
        return 1;
    }
    DftiFreeDescriptor(&hand);
    return 0;
}

/*
 * sizeof(x) == n1*n2*n3*sizeof(double complex)
 */
int ifft3_inplace(double complex *x, int n1, int n2, int n3)
{
    MKL_LONG status = 0;
    DFTI_DESCRIPTOR_HANDLE hand = 0;
    {
        MKL_LONG N[3]; N[0] = n1; N[1] = n2; N[2] = n3;
        status = DftiCreateDescriptor(&hand, DFTI_DOUBLE, DFTI_COMPLEX, 3, N);
        if (0 != status) {
            DftiFreeDescriptor(&hand);
            return 1;
        }
    }
    status = DftiCommitDescriptor(hand);
    if (0 != status) {
        DftiFreeDescriptor(&hand);
        return 1;
    }
    MKL_Complex16 *x_ = x;
    status = DftiComputeBackward(hand, x_);
    if (0 != status) {
        DftiFreeDescriptor(&hand);
        return 1;
    }
    DftiFreeDescriptor(&hand);
    return 0;
}
