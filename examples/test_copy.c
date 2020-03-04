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

 printf("===TEST mncblas_ccopy()===\n");
 printf("c2 avant : [");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("{%.0f %.0f}", c2[idx].real, c2[idx].imaginary);
 }
 printf("]\n");

 mncblas_ccopy(TAILLE, c1, 1, c2, 1);

 printf("c2 apres : [");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("{%.0f %.0f}", c2[idx].real, c2[idx].imaginary);
 }
 printf("]\n");


printf("\n===TEST mncblas_zcopy()===\n");
 printf("cd2 avant : [");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("{%.0f %.0f}", cd2[idx].real, cd2[idx].imaginary);
 }
 printf("]\n");

 mncblas_zcopy(TAILLE, cd1, 1, cd2, 1);

 printf("cd2 apres : [");
 for (int idx=0; idx<TAILLE; idx++) {
     printf("{%.0f %.0f}", cd2[idx].real, cd2[idx].imaginary);
 }
 printf("]\n\n");



 printf("===TEST perf mncblas_scopy===\n");
 float f1[TAILLE];
 float f2[TAILLE];
 init_flop () ;
 start =_rdtsc () ;
 for (i = 0 ; i < NB_FOIS; i++) {
     mncblas_scopy(TAILLE, f1, 1, f2, 1);
 }
 end = _rdtsc () ;
 printf ("apres boucle %lld cycles \n\n", end-start) ;


 printf("===TEST perf mncblas_dcopy===\n");
 double d1[TAILLE];
 double d2[TAILLE];
 init_flop () ;
 start =_rdtsc () ;
 for (i = 0 ; i < NB_FOIS; i++) {
     mncblas_dcopy(TAILLE, d1, 1, d2, 1);
 }
 end = _rdtsc () ;
 printf ("apres boucle %lld cycles \n\n", end-start) ;


 printf("===TEST perf mncblas_ccopy===\n");
 init_flop () ;
 start =_rdtsc () ;
 for (i = 0 ; i < NB_FOIS; i++) {
     mncblas_ccopy(TAILLE, c1, 1, c2, 1);
 }
 end = _rdtsc () ;
 printf ("apres boucle %lld cycles \n\n", end-start) ;


 printf("===TEST perf mncblas_ccopy===\n");
 init_flop () ;
 start =_rdtsc () ;
 for (i = 0 ; i < NB_FOIS; i++) {
     mncblas_ccopy(TAILLE, cd1, 1, cd2, 1);
 }
 end = _rdtsc () ;
 printf ("apres boucle %lld cycles \n\n", end-start) ;

  exit (0) ;
}
