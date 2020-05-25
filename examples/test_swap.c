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
     c1[idx] = (complexe_float_t){(float)idx,(float)idx};
     c2[idx] = (complexe_float_t){0.0,0.0};
     cd1[idx] = (complexe_double_t){(double)idx,(double)idx};
     cd2[idx] = (complexe_double_t){0.0,0.0};
 }

 printf("===TEST mncblas_cswap()===\n");
 printf("c1 avant : [");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("{%.0f %.0f}", c1[idx].real, c1[idx].imaginary);
 }
 printf("]\n");
 printf("c2 avant : [");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("{%.0f %.0f}", c2[idx].real, c2[idx].imaginary);
 }
 printf("]\n");

 mncblas_cswap(TAILLE, c1, 1, c2, 1);

 printf("c1 apres : [");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("{%.0f %.0f}", c1[idx].real, c1[idx].imaginary);
 }
 printf("]\n");
 printf("c2 apres : [");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("{%.0f %.0f}", c2[idx].real, c2[idx].imaginary);
 }
 printf("]\n");


printf("\n===TEST mncblas_zswap()===\n");
 printf("cd1 avant : [");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("{%.0f %.0f}", cd1[idx].real, cd1[idx].imaginary);
 }
 printf("]\n");
 printf("cd2 avant : [");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("{%.0f %.0f}", cd2[idx].real, cd2[idx].imaginary);
 }
 printf("]\n");

 mncblas_zswap(TAILLE, cd1, 1, cd2, 1);

 printf("cd1 apres : [");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("{%.0f %.0f}", cd1[idx].real, cd1[idx].imaginary);
 }
 printf("]\n");
 printf("cd2 apres : [");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("{%.0f %.0f}", cd2[idx].real, cd2[idx].imaginary);
 }
 printf("]\n\n");



 printf("===TEST perf mncblas_sswap===\n");
 float f1[TAILLE];
 float f2[TAILLE];
 init_flop () ;
 start =_rdtsc () ;
 for (i = 0 ; i < NB_FOIS; i++) {
     mncblas_sswap(TAILLE, f1, 1, f2, 1);
 }
 end = _rdtsc () ;
 printf ("apres boucle %lld cycles \n", end-start) ;
 calcul_octet("sswap", TAILLE * NB_FOIS * 3, end-start);
 printf("\n");


 printf("===TEST perf mncblas_dswap===\n");
 double d1[TAILLE];
 double d2[TAILLE];
 init_flop () ;
 start =_rdtsc () ;
 for (i = 0 ; i < NB_FOIS; i++) {
     mncblas_dswap(TAILLE, d1, 1, d2, 1);
 }
 end = _rdtsc () ;
 printf ("apres boucle %lld cycles \n", end-start) ;
 calcul_octet("dswap", TAILLE * NB_FOIS * 3, end-start);
 printf("\n");


 printf("===TEST perf mncblas_cswap===\n");
 init_flop () ;
 start =_rdtsc () ;
 for (i = 0 ; i < NB_FOIS; i++) {
     mncblas_cswap(TAILLE, c1, 1, c2, 1);
 }
 end = _rdtsc () ;
 printf ("apres boucle %lld cycles \n", end-start) ;
 calcul_octet("cswap", TAILLE * NB_FOIS * 6, end-start);
 printf("\n");

 printf("===TEST perf mncblas_zswap===\n");
 init_flop () ;
 start =_rdtsc () ;
 for (i = 0 ; i < NB_FOIS; i++) {
     mncblas_zswap(TAILLE, cd1, 1, cd2, 1);
 }
 end = _rdtsc () ;
 printf ("apres boucle %lld cycles \n", end-start) ;
 calcul_octet("zswap", TAILLE * NB_FOIS * 6, end-start);
 printf("\n");
  exit (0) ;
}
