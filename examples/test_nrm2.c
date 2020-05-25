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

    printf("TEST NRM2\n");

    printf("\n=== mnblas_snrm2 ===\n\n");
    for (i = 0 ; i < NB_FOIS; i++) {
        float res;
        vector_init_s(vef, 1.0);
        res = 0.0;

        start = _rdtsc () ;
            res = mnblas_snrm2(VECSIZE, vef, 1) ;
        end =_rdtsc() ;
        printf ("mnblas_snrm2 %d : res = %3.2f nombre de cycles: %Ld\n", i, res, end-start) ;
        calcul_flop ("snrm2 ", VECSIZE * 2 + 1, end-start) ;
   }

   printf("\n=== mnblas_dnrm2 ===\n\n");
   for (i = 0 ; i < NB_FOIS; i++) {
       double res;
       vector_init_d (ved, 1.0);
       res = 0.0 ;

       start = _rdtsc () ;
           res = mnblas_dnrm2(VECSIZE, ved, 1) ;
       end = _rdtsc () ;

       printf ("mnblas_dnrm2 %d : res = %3.2f nombre de cycles: %Ld\n", i, res, end-start) ;
       calcul_flop ("dnrm2 ", VECSIZE * 2 + 1, end-start) ;
  }

  printf("\n=== mnblas_scnrm2 ===\n\n");
  for (i = 0 ; i < NB_FOIS; i++) {
      float res;
      complexe_float_t val = {1.0,1.0};
      vector_init_c (vecf, val);
      res = 0.0;

      start = _rdtsc () ;
          res = mnblas_scnrm2 (VECSIZE, vecf, 1);
      end = _rdtsc () ;

      printf ("mnblas_scnrm2 %d : res = %3.2f nombre de cycles: %Ld\n", i, res, end-start) ;
      calcul_flop ("scnrm2 ", 1 + 8 * VECSIZE, end-start) ;
 }

 printf("\n=== mnblas_dznrm2 ===\n\n");
 for (i = 0 ; i < NB_FOIS; i++) {
     double res;
     complexe_double_t val = {1.0,1.0};
     vector_init_z (vecd, val) ;
     res = 0.0;

     start = _rdtsc () ;
         res = mnblas_dznrm2 (VECSIZE, vecd, 1) ;
     end = _rdtsc () ;

     printf ("mnblas_dznrm2 %d : res = %3.2f nombre de cycles: %Ld\n", i, res, end-start) ;
     calcul_flop ("dznrm2 ", 1 + 8 * VECSIZE, end-start) ;
 }
}
