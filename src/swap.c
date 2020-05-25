#include "mnblas.h"

void mncblas_sswap(const int N, float *X, const int incX,
                 float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float save ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }

  return ;
}

void mncblas_dswap(const int N, double *X, const int incX,
                 double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register double save ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      save = Y [j] ;
      Y [j] = X [i] ;
      X [i] = save ;
    }
  return ;
}

void mncblas_cswap(const int N, void *X, const int incX,
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float save ;

  float *x = (float *)X;
  float *y = (float *)Y;

  for (; ((i < N*2) && (j < N*2)) ; i += incX, j+=incY)
    {
      save = y[j] ;
      y[j] = x[i] ;
      x[i] = save ;
    }
  return ;
}

void mncblas_zswap(const int N, void *X, const int incX,
		                    void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register double save ;

  double *x = (double *)X;
  double *y = (double *)Y;

  for (; ((i < N*2) && (j < N*2)) ; i += incX, j+=incY)
    {
      save = y[j] ;
      y[j] = x[i] ;
      x[i] = save ;
    }
  return ;
}
