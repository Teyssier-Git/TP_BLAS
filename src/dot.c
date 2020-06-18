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
  float dot = 0.0 ;

  int min;
  if (N/incX + 1 < N/incY + 1) {
      min = N/incX+1;
  } else {
      min = N/incY+1;
  }

  #pragma omp parallel for reduction(+:dot)
  for (i = 0 ; i < min ; i ++)
    {
      dot += X [i*incX] * Y [i*incY] ;
    }

  return dot;
}

double mncblas_ddot(const int N, const double *X, const int incX,
                 const double *Y, const int incY)
{
    register unsigned int i = 0 ;
    double dot = 0.0 ;

    int min;
    if (N/incX + 1 < N/incY + 1) {
        min = N/incX+1;
    } else {
        min = N/incY+1;
    }

    #pragma omp parallel for reduction(+:dot)
    for (i = 0 ; i < min ; i++)
      {
        dot += X [i*incX] * Y [i*incY] ;
      }

    return dot ;
}

void   mncblas_cdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
    register unsigned int i = 0 ;

    int min;
    if (N/incX + 1 < N/incY + 1) {
        min = N/incX+1;
    } else {
        min = N/incY+1;
    }

    complexe_float_t *fX = (complexe_float_t *)X;
    complexe_float_t *fY = (complexe_float_t *)Y;
    complexe_float_t dot;
    dot.real = 0.0;
    dot.imaginary = 0.0;

    #pragma omp declare reduction(add_cf : complexe_float_t : omp_out = add_complexe_float(omp_out, omp_in))
    #pragma omp parallel for reduction(add_cf:dot)
    for (i = 0 ; i < min ; i++)
      {
        dot = add_complexe_float(dot, mult_complexe_float(fX[i*incX], fY[i*incY]));
      }

      complexe_float_t *rep = (complexe_float_t *)dotu;
      rep->real = dot.real;
      rep->imaginary = dot.imaginary;
}

void   mncblas_cdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
    register unsigned int i = 0 ;

    int min;
    if (N/incX + 1 < N/incY + 1) {
        min = N/incX+1;
    } else {
        min = N/incY+1;
    }

    complexe_float_t *fX = (complexe_float_t *)X;
    complexe_float_t *fY = (complexe_float_t *)Y;
    complexe_float_t dot;
    dot.real = 0.0;
    dot.imaginary = 0.0;

    #pragma omp declare reduction(add_cf : complexe_float_t : omp_out = add_complexe_float(omp_out, omp_in))
    #pragma omp parallel for reduction(add_cf:dot)
    for (i = 0 ; i < min ; i++)
      {
        dot = add_complexe_float(dot, mult_complexe_float(conj_complexe_float(fX[i*incX]), conj_complexe_float(fY[i*incY])));
      }

      complexe_float_t *rep = (complexe_float_t *)dotc;
      rep->real = dot.real;
      rep->imaginary = dot.imaginary;
}

void   mncblas_zdotu_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotu)
{
    register unsigned int i = 0 ;

    int min;
    if (N/incX + 1 < N/incY + 1) {
        min = N/incX+1;
    } else {
        min = N/incY+1;
    }

    complexe_double_t *fX = (complexe_double_t *)X;
    complexe_double_t *fY = (complexe_double_t *)Y;
    complexe_double_t dot;
    dot.real = 0.0;
    dot.imaginary = 0.0;

    #pragma omp declare reduction(add_cd : complexe_double_t : omp_out = add_complexe_double(omp_out, omp_in))
    #pragma omp parallel for reduction(add_cd:dot)
    for (i = 0 ; i < min ; i++)
      {
        dot = add_complexe_double(dot, mult_complexe_double(fX[i*incX], fY[i*incY]));
      }

      complexe_double_t *rep = (complexe_double_t *)dotu;
      rep->real = dot.real;
      rep->imaginary = dot.imaginary;
}

void   mncblas_zdotc_sub(const int N, const void *X, const int incX,
                       const void *Y, const int incY, void *dotc)
{
    register unsigned int i = 0 ;

    int min;
    if (N/incX + 1 < N/incY + 1) {
        min = N/incX+1;
    } else {
        min = N/incY+1;
    }

    complexe_double_t *fX = (complexe_double_t *)X;
    complexe_double_t *fY = (complexe_double_t *)Y;
    complexe_double_t dot;
    dot.real = 0.0;
    dot.imaginary = 0.0;

    #pragma omp declare reduction(add_cd : complexe_double_t : omp_out = add_complexe_double(omp_out, omp_in))
    #pragma omp parallel for reduction(add_cd:dot)
    for (i = 0 ; i < min ; i++) {
        dot = add_complexe_double(dot, mult_complexe_double(conj_complexe_double(fX[i*incX]), conj_complexe_double(fY[i*incY])));
    }

    complexe_double_t *rep = (complexe_double_t *)dotc;
    rep->real = dot.real;
    rep->imaginary = dot.imaginary;
}
