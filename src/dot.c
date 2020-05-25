#include "mnblas.h"
#include "../include/complexe.h"
#include <stdio.h>

/*
float mncblas_sdot(const int N, const float *X, const int incX,
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  register float dot = 0.0 ;

  for (; ((i < N) && (j < N)) ; i += incX, j+=incY)
    {
      dot = dot + X [i] * Y [j] ;
    }

  return dot ;
}
*/

float mncblas_sdot(const int N, const float *X, const int incX,
                 const float *Y, const int incY)
{
  register unsigned int i = 0 ;
  register unsigned int j = 0 ;
  float dot = 0.0 ;


  for (i = 0 ; i < N ; i += incX)
    {
      dot += X [i] * Y [j] ;
      j+=incY ;
    }

  return dot ;
}

double mncblas_ddot(const int N, const double *X, const int incX,
                 const double *Y, const int incY)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;
    double dot = 0.0 ;


    for (i = 0 ; i < N ; i += incX)
      {
        dot += X [i] * Y [j] ;
        j+=incY ;
      }

    return dot ;
}

void   mncblas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;

    complexe_float_t *fX = (complexe_float_t *)X;
    complexe_float_t *fY = (complexe_float_t *)Y;
    complexe_float_t *dot = (complexe_float_t *)dotu;
    dot->real = 0.0;
    dot->imaginary = 0.0;

    for (i = 0 ; i < N ; i += incX)
      {
        *dot = add_complexe_float(*dot, mult_complexe_float(fX[i], fY[j]));
        j+=incY ;
      }
}

void   mncblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;

    complexe_float_t *fX = (complexe_float_t *)X;
    complexe_float_t *fY = (complexe_float_t *)Y;
    complexe_float_t *dot = (complexe_float_t *)dotc;
    dot->real = 0.0;
    dot->imaginary = 0.0;

    for (i = 0 ; i < N ; i += incX)
      {
        *dot = add_complexe_float(*dot, mult_complexe_float(conj_complexe_float(fX[i]), conj_complexe_float(fY[j])));
        j+=incY ;
      }
}

void   mncblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;

    complexe_double_t *fX = (complexe_double_t *)X;
    complexe_double_t *fY = (complexe_double_t *)Y;
    complexe_double_t *dot = (complexe_double_t *)dotu;
    dot->real = 0.0;
    dot->imaginary = 0.0;

    for (i = 0 ; i < N ; i += incX)
      {
        *dot = add_complexe_double(*dot, mult_complexe_double(fX[i], fY[j]));
        j+=incY ;
      }
}

void   mncblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
    register unsigned int i = 0 ;
    register unsigned int j = 0 ;

    complexe_double_t *fX = (complexe_double_t *)X;
    complexe_double_t *fY = (complexe_double_t *)Y;
    complexe_double_t *dot = (complexe_double_t *)dotc;
    dot->real = 0.0;
    dot->imaginary = 0.0;

    for (i = 0 ; i < N ; i += incX)
      {
        *dot = add_complexe_double(*dot, mult_complexe_double(conj_complexe_double(fX[i]), conj_complexe_double(fY[j])));
        j+=incY ;
      }
}
