#include "mnblas.h"

void mncblas_scopy(const int N, const float *X, const int incX,
                 float *Y, const int incY)
{
    register unsigned int i;

    int min;
    if (N/incX + 1 < N/incY + 1) {
        min = N/incX+1;
    } else {
        min = N/incY+1;
    }

    #pragma omp parallel for
    for (i=0; i < min; i++) {
      Y [i*incY] = X [i*incX] ;
    }
}

void mncblas_dcopy(const int N, const double *X, const int incX,
                 double *Y, const int incY)
{
    register unsigned int i;

    int min;
    if (N/incX + 1 < N/incY + 1) {
        min = N/incX+1;
    } else {
        min = N/incY+1;
    }

    #pragma omp parallel for
    for (i=0; i < min; i++) {
      Y [i*incY] = X [i*incX] ;
    }
}

void mncblas_ccopy(const int N, const void *X, const int incX,
		                    void *Y, const int incY)
{
    float *y = (float *)Y;
    float *x = (float *)X;

    register unsigned int i;

    int min;
    if (2*N/incX + 1 < 2*N/incY + 1) {
        min = 2*N/incX+1;
    } else {
        min = 2*N/incY+1;
    }

    #pragma omp parallel for
    for (i=0; i < min; i++) {
      y[i*incY] = x[i*incX] ;
    }
}

void mncblas_zcopy(const int N, const void *X, const int incX,
		                    void *Y, const int incY)
{
    double *y = (double *)Y;
    double *x = (double *)X;

    register unsigned int i;

    int min;
    if (2*N/incX + 1 < 2*N/incY + 1) {
        min = 2*N/incX+1;
    } else {
        min = 2*N/incY+1;
    }

    #pragma omp parallel for
    for (i=0; i < min; i++) {
      y[i*incY] = x[i*incX] ;
    }
}
