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

vfloat vef1, vef2;
vdouble ved1, ved2;
vCompFloat vecf1, vecf2;
vCompDouble vecd1, vecd2;

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

    printf("TEST DOT\n");

    printf("\n=== mncblas_sdot ===\n\n");
    for (i = 0 ; i < NB_FOIS; i++) {
        float res;
        vector_init_s (vef1, 1.0) ;
        vector_init_s (vef2, 2.0) ;
        res = 0.0 ;

        start = _rdtsc () ;
            res = mncblas_sdot (VECSIZE, vef1, 1, vef2, 1) ;
        end = _rdtsc () ;
//      "\033[31;01m%c\033[00m"
        printf ("\033[33;01mmncblas_sdot %d :\033[00m \033[35;01mres = %3.2f nombre de cycles: %Ld\033[00m\n", i, res, end-start) ;
        calcul_flop ("sdot ", 2 * VECSIZE, end-start) ;
   }

   printf("\n=== mncblas_ddot ===\n\n");
   for (i = 0 ; i < NB_FOIS; i++) {
       double res;
       vector_init_d (ved1, 1.0) ;
       vector_init_d (ved2, 2.0) ;
       res = 0.0 ;

       start = _rdtsc () ;
           res = mncblas_ddot (VECSIZE, ved1, 1, ved2, 1) ;
       end = _rdtsc () ;

       printf ("\033[33;01mmncblas_sdot %d :\033[00m \033[35;01mres = %3.2f nombre de cycles: %Ld\033[00m\n", i, res, end-start) ;
       calcul_flop ("ddot ", 2 * VECSIZE, end-start) ;
  }

  printf("\n=== mncblas_cdotu_sub ===\n\n");
  for (i = 0 ; i < NB_FOIS; i++) {
      complexe_float_t res;
      complexe_float_t val1 = {1.0,1.0};
      complexe_float_t val2 = {2.0,2.0};
      vector_init_c (vecf1, val1) ;
      vector_init_c (vecf2, val2) ;
      res.real = 0.0 ;
      res.imaginary = 0.0 ;

      start = _rdtsc () ;
          mncblas_cdotu_sub (VECSIZE, vecf1, 1, vecf2, 1, &res) ;
      end = _rdtsc () ;

      printf ("\033[33;01mmncblas_cdotu_sub %d :\033[00m \033[35;01mres = %3.2f + i*%3.2f nombre de cycles: %Ld\033[00m\n", i, res.real, res.imaginary, end-start) ;
      calcul_flop ("cdotu_sub ", 8 * VECSIZE, end-start) ;
 }

 printf("\n=== mncblas_zdotu_sub ===\n\n");
 for (i = 0 ; i < NB_FOIS; i++) {
     complexe_double_t res;
     complexe_double_t val1 = {1.0,1.0};
     complexe_double_t val2 = {2.0,2.0};
     vector_init_z (vecd1, val1) ;
     vector_init_z (vecd2, val2) ;
     res.real = 0.0 ;
     res.imaginary = 0.0 ;

     start = _rdtsc () ;
         mncblas_zdotu_sub (VECSIZE, vecd1, 1, vecd2, 1, &res) ;
     end = _rdtsc () ;

     printf ("\033[33;01mmncblas_zdotu_sub %d :\033[00m \033[35;01mres = %3.2f + i*%3.2f nombre de cycles: %Ld\033[00m\n", i, res.real, res.imaginary, end-start) ;
     calcul_flop ("zdotu_sub ", 8 * VECSIZE, end-start) ;
 }
}
