#include <stdio.h>
#include <x86intrin.h>

#include "../include/mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE    10

#define NB_FOIS    10

typedef float vfloat [VECSIZE] ;
typedef complexe_double_t vCompDouble [VECSIZE];

vfloat vec1, vec2 ;
vCompDouble vecd1, vecd2, valpha;

void vector_init (vfloat V, float x)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    V [i] = x ;

  return ;
}

void vector_init_z (vCompDouble V, complexe_double_t x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE; i++) {
        V[i].real = x.real;
        V[i].imaginary = x.imaginary;
    }
    return;
}


void vector_print (vfloat V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f ", V[i]) ;
  printf ("\n") ;

  return ;
}

void vector_print_z (vCompDouble V)
{
  register unsigned int i ;

  for (i = 0; i < VECSIZE; i++)
    printf ("%f+i%f ", V[i].real, V[i].imaginary) ;
  printf ("\n") ;

  return ;
}


int main (int argc, char **argv)
{
 unsigned long long start, end ;
 int i ;

 init_flop () ;

 for (i = 0 ; i < NB_FOIS; i++)
   {
     vector_init (vec1, 2.0) ;
     vector_init (vec2, 4.0) ;
     start = _rdtsc () ;
        mnblas_saxpy(VECSIZE, 3.0, vec1,1, vec2, 1) ;
     end = _rdtsc () ;

     printf ("mnblas_saxpy %d : nombre de cycles: %Ld \n", i, end-start) ;
     vector_print(vec2);
     calcul_flop ("saxpy ", 2 * VECSIZE, end-start) ;
   }
  for (i = 0 ; i < NB_FOIS; i++) {

    complexe_double_t val1 = {2.0,2.0};
    complexe_double_t val2 = {4.0,4.0};
    complexe_double_t alpha = {3.0, 0.0};
    vector_init_z (vecd1, val1) ;
    vector_init_z (vecd2, val2) ;
    vector_init_z (valpha, alpha) ;


    start = _rdtsc () ;
        mnblas_zaxpy (VECSIZE, valpha, vecd1, 1, vecd2, 1) ;
    end = _rdtsc () ;

    printf ("mnblas_zaxpy %d :  nombre de cycles: %Ld\n", i, end-start) ;
    vector_print_z(vecd2);
    calcul_flop ("zaxpy ", 10 * VECSIZE, end-start) ;
}

}
