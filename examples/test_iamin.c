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

vfloat vef1;
vdouble ved1;
vCompFloat vecf1;
vCompDouble vecd1;

void vector_init_s (vfloat V, float x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE; i++) {
        V [i] = x + (12*i) % 21 ;
    }
    return;
}

void vector_init_d (vdouble V, double x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE; i++) {
       V [i] = x + (12*i) % 21 ;
    }
    return;
}

void vector_init_c (vCompFloat V, complexe_float_t x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE; i++) {
      V [i].real = x.real + (12*i) % 21 ;
      V [i].imaginary = x.imaginary + (12*i) % 22 ;
    }
    return;
}

void vector_init_z (vCompDouble V, complexe_double_t x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE; i++) {
      V [i].real = x.real + (12*i) % 21 ;
      V [i].imaginary = x.imaginary + (12*i) % 22 ;
    }
    return;
}

void vfloat_print (vfloat V) {
    register unsigned int i;
    for (i = 0; i < VECSIZE; i++) {
        printf ("%d : %f ", i, V[i]) ;
    }
    printf ("\n") ;
    return;
}

void vdouble_print (vdouble V){
    register unsigned int i;
    for (i = 0; i < VECSIZE; i++) {
        printf ("%d : %lf ", i, V[i]) ;
    }
    printf ("\n") ;
    return;
}

void vCompFloat_print (vCompFloat V){
    register unsigned int i;
    for (i = 0; i < VECSIZE; i++) {
        printf ("%d : %f + i %f", i, V[i].real, V[i].imaginary) ;
    }
    printf ("\n") ;
    return;
}

void vCompDouble_print (vCompDouble V){
    register unsigned int i;
    for (i = 0; i < VECSIZE; i++) {
        printf ("%d : %lf + i %lf |", i, V[i].real, V[i].imaginary) ;
    }
    printf ("\n") ;
    return;
}

int main (int argc, char **argv) {
    unsigned long long start, end;
    int i;

    init_flop () ;

    printf("TEST DOT\n");

    printf("\n=== mnblas_isamin ===\n\n");
    for (i = 0 ; i < NB_FOIS; i++) {
        int res;
        vector_init_s (vef1, 1.0) ;
        vfloat_print(vef1);

        start = _rdtsc () ;
            res = mnblas_isamin (VECSIZE, vef1, 1) ;
        end = _rdtsc () ;
        
        printf ("mnblas_isamin %d : res = %d nombre de cycles: %Ld\n", i, res, end-start) ;
        calcul_flop ("isamin ", 2 * VECSIZE, end-start) ;
   }

   printf("\n=== mnblas_idamin ===\n\n");
   for (i = 0 ; i < NB_FOIS; i++) {
       int res;
       vector_init_d (ved1, 1.0) ;
       vdouble_print(ved1);

       start = _rdtsc () ;
           res = mnblas_idamin (VECSIZE, ved1, 1) ;
       end = _rdtsc () ;

       printf ("mnblas_idamin %d : res = %d nombre de cycles: %Ld\n", i, res, end-start) ;
       calcul_flop ("idamin ", 2 * VECSIZE, end-start) ;
  }

  printf("\n=== mnblas_icamin ===\n\n");
  for (i = 0 ; i < NB_FOIS; i++) {
      int res;
      complexe_float_t val1 = {1.0,1.0};
      vector_init_c (vecf1, val1) ;
      vCompFloat_print(vecf1);

      start = _rdtsc () ;
        res = mnblas_icamin (VECSIZE, vecf1, 1) ;
      end = _rdtsc () ;

      printf ("mnblas_icamin %d : res = %d nombre de cycles: %Ld\n", i, res, end-start) ;
      calcul_flop ("icamin ", 8 * VECSIZE, end-start) ;
 }

 printf("\n=== mnblas_izamin ===\n\n");
 for (i = 0 ; i < NB_FOIS; i++) {
     int res;
     complexe_double_t val1 = {1.0,1.0};
     vector_init_z (vecd1, val1) ;
     vCompDouble_print(vecd1);

     start = _rdtsc () ;
         res = mnblas_izamin (VECSIZE, vecf1, 1) ;
     end = _rdtsc () ;

     printf ("mnblas_izamin %d : res = %d nombre de cycles: %Ld\n", i, res, end-start) ;
     calcul_flop ("izamin ", 10 * VECSIZE, end-start) ;
}


}
