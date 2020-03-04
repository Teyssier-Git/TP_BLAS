#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define    NB_FOIS        4194304
#define    TAILLE         10

int main (int argc, char **argv)
{
 complexe_float_t c1[TAILLE];
 complexe_float_t c2[TAILLE];

 // complexe_double_t cd1[TAILLE] ;
 // complexe_double_t cd2[TAILLE] ;

 // unsigned long long int start, end ;
 // int i ;

 for (int idx = 0; idx<TAILLE; idx++) {
     c1[idx] = (complexe_float_t){0.0,0.0};
 }

 printf("c2\n");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("%.0f %.0f\n", c2[idx].real, c2[idx].imaginary);
 }
 mncblas_ccopy(TAILLE, c1, 1, c2, 1);
 printf("\nc2\n");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("%.0f %.0f\n", c2[idx].real, c2[idx].imaginary);
 }
 // init_flop () ;
 //
 // printf ("cd3.r %f cd3.i %f\n", cd3.real, cd3.imaginary) ;
 //
 // start =_rdtsc () ;
 //
 // for (i = 0 ; i < NB_FOIS; i++)
 //   {
 //     mncblas_ccopy(const int N, const void *X, const int incX,
 //     		                    void *Y, const int incY)
 //   }
 //
 // end = _rdtsc () ;
 //
 //  printf ("apres boucle cd3.real %f cd3.imaginary %f %lld cycles \n", cd3.real, cd3.imaginary, end-start) ;
 //
 //  calcul_flop ("calcul copy ", NB_FOIS*4, end-start) ;
  exit (0) ;
}
