#include "../include/mnblas.h"
#include "../include/complexe.h"
#include <stdio.h>

void mnblas_saxpy(const int N, const float alpha, const float *X,
                 const int incX, float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;


  for ( i=0; i < N ; i += incX )
    {
      Y[j] = X [i] * alpha + Y[j] ;
      j+=incY ;
    }

  return;
}

void mnblas_daxpy(const int N, const double alpha, const double *X,
                 const int incX, double *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;

  for ( i=0; i < N ; i += incX )
    {
      Y[j] = X [i] * alpha + Y[j] ;
      j+=incY ;
    }

  return;
}

void mnblas_caxpy(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_float_t *fX = (complexe_float_t *)X;
  complexe_float_t *fY = (complexe_float_t *)Y;
  complexe_float_t *falpha = (complexe_float_t *)alpha;

  for ( i=0; i < N ; i += incX )
    {
      fY[j]=add_complexe_float(fY[j],mult_complexe_float(fX[i],falpha[i])) ;
      j+=incY ;
    }

  return;
}

void mnblas_zaxpy(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  complexe_double_t *dX = (complexe_double_t *)X;
  complexe_double_t *dY = (complexe_double_t *)Y;
  complexe_double_t *dalpha = (complexe_double_t *)alpha;

  for ( i=0; i < N ; i += incX )
    {
      dY[j]=add_complexe_double(dY[j],mult_complexe_double(dX[i],dalpha[i])) ;
      j+=incY ;
    }

  return;
}
