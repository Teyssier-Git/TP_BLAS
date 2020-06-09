#include "mnblas.h"
#include "complexe.h"
#include <stdio.h>
#include <math.h>

float  mnblas_snrm2(const int N, const float *X, const int incX) {
    float res = 0.0;

    register unsigned int i;

    for (i = 0 ; i < N ; i += incX) {
        res += X[i]*X[i];
    }

    return sqrt(res);
}

double mnblas_dnrm2(const int N, const double *X, const int incX) {
    double res = 0.0;

    register unsigned int i;

    for (i = 0 ; i < N ; i += incX) {
        res += X[i]*X[i];
    }

    return sqrt(res);
}

float  mnblas_scnrm2(const int N, const void *X, const int incX) {
    float res = 0.0;

    register unsigned int i;
    complexe_float_t *x = (complexe_float_t *)X;

    for (i = 0 ; i < N ; i += incX) {
        complexe_float_t temp = mult_complexe_float(x[i], conj_complexe_float(x[i]));
        res += temp.real;
    }

    return sqrt(res);
}

double mnblas_dznrm2(const int N, const void *X, const int incX) {
    double res = 0.0;

    register unsigned int i;
    complexe_double_t *x = (complexe_double_t *)X;

    for (i = 0 ; i < N ; i += incX) {
        complexe_double_t temp = mult_complexe_double(x[i], conj_complexe_double(x[i]));
        res += temp.real;
    }

    return sqrt(res);
}
