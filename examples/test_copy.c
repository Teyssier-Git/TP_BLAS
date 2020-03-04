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

 complexe_double_t cd1[TAILLE] ;
 complexe_double_t cd2[TAILLE] ;

 unsigned long long int start, end ;
 int i ;

 for (int idx = 0; idx<TAILLE; idx++) {
     c1[idx] = (complexe_float_t){0.0,0.0};
     cd1[idx] = (complexe_double_t){0.0,0.0};
     c2[idx] = (complexe_float_t){1.0,1.0};
     cd2[idx] = (complexe_double_t){1.0,1.0};
 }

 for (int idx = 0; idx<TAILLE; idx++) {

 }

 printf("c2 avant\n [");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("{%.0f %.0f}", c2[idx].real, c2[idx].imaginary);
 }
 printf("]");

 mncblas_ccopy(TAILLE, c1, 1, c2, 1);

 printf("\nc2 apres\n [");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("{%.0f %.0f}", c2[idx].real, c2[idx].imaginary);
 }
 printf("]");

 printf("cd2 avant\n [");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("{%.0f %.0f}", cd2[idx].real, cd2[idx].imaginary);
 }
 printf("]");

 mncblas_zcopy(TAILLE, cd1, 1, cd2, 1);

 printf("\ncd2 apres\n [");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("{%.0f %.0f}", c2[idx].real, c2[idx].imaginary);
 }
 printf("]\n\n\n");

 init_flop () ;

 start =_rdtsc () ;

 for (i = 0 ; i < NB_FOIS; i++)
   {
     mncblas_ccopy(TAILLE, c1, 1, c2, 1);
   }

 end = _rdtsc () ;

  printf ("apres boucle cd3.real %f cd3.imaginary %f %lld cycles \n", c1[0].real, c1[0].imaginary, end-start) ;

  calcul_flop ("calcul copy ", NB_FOIS*4, end-start) ;
  exit (0) ;
}
