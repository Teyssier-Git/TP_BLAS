#include "mnblas.h"
#include "complexe.h"
#include <stdio.h>

float  mnblas_sasum(const int N, const float *X, const int incX) {
    float res = 0.0;

    register unsigned int i;

    for (i = 0 ; i < N ; i += incX) {
        if (X[i] >= 0) {
            res += X [i];
        } else {
            res -= X [i];
        }
    }

    return res;
}

double mnblas_dasum(const int N, const double *X, const int incX) {
    double res = 0.0;

    register unsigned int i;

    for (i = 0 ; i < N ; i += incX) {
        if (X[i] >= 0) {
            res += X [i];
        } else {
            res -= X [i];
        }
    }

    return res;
}

float  mnblas_scasum(const int N, const void *X, const int incX) {
    float res = 0.0;

    register unsigned int i;
    float *x = (float *)X;

    for (i = 0 ; i < N*2 ; i += incX) {
        if (x[i] >= 0) {
            res += x[i];
        } else {
            res -= x[i];
        }
    }

    return res;
}

double mnblas_dzasum(const int N, const void *X, const int incX) {
    double res = 0.0;

    register unsigned int i;
    double *x = (double *)X;

    for (i = 0 ; i < N*2 ; i += incX) {
        if (x[i] >= 0) {
            res += x[i];
        } else {
            res -= x[i];
        }
    }

    return res;
}
