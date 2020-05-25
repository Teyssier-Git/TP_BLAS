#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE    65536

#define NB_FOIS    10

typedef float vfloat [VECSIZE];
typedef double vdouble [VECSIZE];
typedef complexe_float_t vCompFloat [VECSIZE];
typedef complexe_double_t vCompDouble [VECSIZE];

vfloat vef;
vdouble ved;
vCompFloat vecf;
vCompDouble vecd;

void vector_init_s (vfloat V, float x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE; i++) {
        V [i] = x ;
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

void vector_init_c (vCompFloat V, complexe_float_t x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE; i++) {
        V[i].real = x.real;
        V[i].imaginary = x.imaginary;
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

    printf("TEST ASUM\n");

    printf("\n=== mnblas_sasum ===\n\n");
    for (i = 0 ; i < NB_FOIS; i++) {
        float res;
        vector_init_s(vef, 1.0);
        res = 0.0;

        start = _rdtsc () ;
            res = mnblas_sasum(VECSIZE, vef, 1) ;
        end =_rdtsc() ;
        printf ("mnblas_sasum %d : res = %3.2f nombre de cycles: %Ld\n", i, res, end-start) ;
        calcul_flop ("sasum ", VECSIZE, end-start) ;
   }

   printf("\n=== mnblas_dasum ===\n\n");
   for (i = 0 ; i < NB_FOIS; i++) {
       double res;
       vector_init_d (ved, 1.0);
       res = 0.0 ;

       start = _rdtsc () ;
           res = mnblas_dasum(VECSIZE, ved, 1) ;
       end = _rdtsc () ;

       printf ("mnblas_dasum %d : res = %3.2f nombre de cycles: %Ld\n", i, res, end-start) ;
       calcul_flop ("dasum ", VECSIZE, end-start) ;
  }

  printf("\n=== mnblas_scasum ===\n\n");
  for (i = 0 ; i < NB_FOIS; i++) {
      float res;
      complexe_float_t val = {1.0,1.0};
      vector_init_c (vecf, val);
      res = 0.0;

      start = _rdtsc () ;
          res = mnblas_scasum (VECSIZE, vecf, 1);
      end = _rdtsc () ;

      printf ("mnblas_cdotu_sub %d : res = %3.2f nombre de cycles: %Ld\n", i, res, end-start) ;
      calcul_flop ("cdotu_sub ", 2 * VECSIZE, end-start) ;
 }

 printf("\n=== mnblas_dzasum ===\n\n");
 for (i = 0 ; i < NB_FOIS; i++) {
     double res;
     complexe_double_t val = {1.0,1.0};
     vector_init_z (vecd, val) ;
     res = 0.0;

     start = _rdtsc () ;
         res = mnblas_dzasum (VECSIZE, vecd, 1) ;
     end = _rdtsc () ;

     printf ("mnblas_zdotc_sub %d : res = %3.2f nombre de cycles: %Ld\n", i, res, end-start) ;
     calcul_flop ("dzasum ", 2 * VECSIZE, end-start) ;
 }
}
