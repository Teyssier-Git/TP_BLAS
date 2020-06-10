#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE    6553

#define NB_FOIS    10

typedef float vfloat [VECSIZE];
typedef float mfloat [VECSIZE*VECSIZE];
typedef double vdouble [VECSIZE];
typedef double mdouble [VECSIZE*VECSIZE];
typedef complexe_float_t vCompFloat [VECSIZE];
typedef complexe_float_t mCompFloat [VECSIZE*VECSIZE];
typedef complexe_double_t vCompDouble [VECSIZE];
typedef complexe_double_t mCompDouble [VECSIZE*VECSIZE];

vfloat vef;
mfloat mef;
vdouble ved;
mdouble med;
vCompFloat vecf;
mCompFloat mecf;
vCompDouble vecd;
mCompDouble mecd;

void vector_init_s (vfloat V, float x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE; i++) {
        V [i] = x ;
    }
    return;
}

void matrix_init_s (mfloat M, float x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE*VECSIZE; i++) {
        M [i] = x ;
    }
    return;
}



void vector_init_d (vdouble V, double x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE; i++) {
        V [i] = x ;
    }
    return;
}

void matrix_init_d (mdouble M, double x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE*VECSIZE; i++) {
        M [i] = x ;
    }
    return;
}



void vector_init_c (vCompFloat V, complexe_float_t x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE; i++) {
        V[i].real = x.real;
        V[i].imaginary = x.imaginary;
    }
    return;
}

void matrix_init_c (mCompFloat M, complexe_float_t x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE*VECSIZE; i++) {
        M[i].real = x.real;
        M[i].imaginary = x.imaginary;
    }
    return;
}



void vector_init_z (vCompDouble V, complexe_double_t x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE; i++) {
        V[i].real = x.real;
        V[i].imaginary = x.imaginary;
    }
    return;
}

void matrix_init_z (mCompDouble M, complexe_double_t x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE*VECSIZE; i++) {
        M[i].real = x.real;
        M[i].imaginary = x.imaginary;
    }
    return;
}



void vector_print (vfloat V) {
    register unsigned int i;
    for (i = 0; i < VECSIZE; i++) {
        printf ("%f ", V[i]) ;
    }
    printf ("\n") ;
    return;
}

int main (int argc, char **argv) {
    unsigned long long start, end;
    int i;

    init_flop () ;

    printf("TEST GEMV\n");

    printf("\n=== mncblas_sgemv ===\n\n");
    for (i = 0 ; i < NB_FOIS; i++) {
        start = _rdtsc () ;
            mncblas_sgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, VECSIZE, 1.0, mef, VECSIZE, vef, 1, 1.0, vef, 1);
        end =_rdtsc() ;
        printf ("mncblas_sgemv %d : nombre de cycles: %Ld\n", i, end-start) ;
        calcul_flop ("sgemv ", VECSIZE * (3 + VECSIZE * 2), end-start) ;
   }

   printf("\n=== mncblas_dgemv ===\n\n");
   for (i = 0 ; i < NB_FOIS; i++) {
       start = _rdtsc () ;
           mncblas_dgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, VECSIZE, 1.0, med, VECSIZE, ved, 1, 1.0, ved, 1) ;
       end = _rdtsc () ;

       printf ("mncblas_dgemv %d : nombre de cycles: %Ld\n", i, end-start) ;
       calcul_flop ("dgemv ", VECSIZE * (3 + VECSIZE * 2), end-start) ;
  }

  printf("\n=== mncblas_cgemv ===\n\n");
  complexe_float_t scal = (complexe_float_t){1.0,0.0};
  for (i = 0 ; i < NB_FOIS; i++) {
      start = _rdtsc () ;
          mncblas_cgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, VECSIZE, &scal, mecf, VECSIZE, vecf, 1, &scal, vecf, 1);
      end = _rdtsc () ;

      printf ("mncblas_cgemv %d : nombre de cycles: %Ld\n", i, end-start) ;
      calcul_flop ("cgemv ", VECSIZE * (14 + VECSIZE * 8), end-start) ;
 }

 printf("\n=== mncblas_zgemv ===\n\n");
 complexe_double_t sal = (complexe_double_t){1.0,0.0};
 for (i = 0 ; i < NB_FOIS; i++) {
     start = _rdtsc () ;
         mncblas_zgemv (MNCblasRowMajor, MNCblasNoTrans, VECSIZE, VECSIZE, &sal, mecd, VECSIZE, vecd, 1, &sal, vecd, 1) ;
     end = _rdtsc () ;

     printf ("mncblas_zgemv %d : nombre de cycles: %Ld\n", i, end-start) ;
     calcul_flop ("zgemv ", VECSIZE * (14 + VECSIZE * 8), end-start) ;
 }
}
