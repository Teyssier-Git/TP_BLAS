#include "mnblas.h"

void mncblas_scopy(const int N, const float *X, const int incX,
                 float *Y, const int incY)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;

    for (; ((i < N) && (j < N)) ; i += incX, j += incY)
    {
      Y [j] = X [i] ;
    }

    return ;
}

void mncblas_dcopy(const int N, const double *X, const int incX,
                 double *Y, const int incY)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;

    for (; ((i < N) && (j < N)) ; i += incX, j += incY)
      {
        Y [j] = X [i] ;
      }

    return ;
}

void mncblas_ccopy(const int N, const void *X, const int incX,
		                    void *Y, const int incY)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;

    float *y = (float *)Y;
    float *x = (float *)X;

    for (; ((i < N*2) && (j < N*2)) ; i += incX, j += incY)
      {
        y [j] = x [i] ;
      }

    return ;
}

void mncblas_zcopy(const int N, const void *X, const int incX,
		                    void *Y, const int incY)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;

    double *y = (double *)Y;
    double *x = (double *)X;

    for (; ((i < N*2) && (j < N*2)) ; i += incX, j += incY)
      {
        y [j] = x [i] ;
      }

    return ;
}
