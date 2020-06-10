#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE    200

#define NB_FOIS    10

typedef float mfloat [VECSIZE*VECSIZE];
typedef double mdouble [VECSIZE*VECSIZE];
typedef complexe_float_t mCompFloat [VECSIZE*VECSIZE];
typedef complexe_double_t mCompDouble [VECSIZE*VECSIZE];

mfloat mef;
mfloat mef2;
mfloat mef3;
mdouble med;
mdouble med2;
mdouble med3;
mCompFloat mecf;
mCompFloat mecf2;
mCompFloat mecf3;
mCompDouble mecd;
mCompDouble mecd2;
mCompDouble mecd3;



void matrix_init_s (mfloat M, float x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE*VECSIZE; i++) {
        M [i] = x ;
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



void matrix_init_c (mCompFloat M, complexe_float_t x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE*VECSIZE; i++) {
        M[i].real = x.real;
        M[i].imaginary = x.imaginary;
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




int main (int argc, char **argv) {
    unsigned long long start, end;
    int i;

    init_flop () ;

    printf("TEST GEMM\n");

       printf("\n=== mncblas_sgemm ===\n\n");
       if (1) {
           float matA[4*4], matB[4*4], matC[4*4];
           for (int i=0; i<4; i++) {
               for (int j=0; j<4; j++) {
                   matA[4*i+j] = i+1;
                   matB[4*i+j] = i+1;
                   matC[4*i+j] = i+1;
               }
           }
           mncblas_sgemm(MNCblasRowMajor, MNCblasNoTrans,MNCblasNoTrans, 4, 4,4, 1.0, matA,4, matB, 4,1.0, matC, 4);

           printf("      | %.1f %.1f %.1f %.1f |   \n",  matC[0], matC[1], matC[2], matC[3]);
           printf("      | %.1f %.1f %.1f %.1f |   \n",  matC[4], matC[5], matC[6], matC[7]);
           printf("      | %.1f %.1f %.1f %.1f |   \n",  matC[8], matC[9], matC[10], matC[11]);
           printf("      | %.1f %.1f %.1f %.1f |   \n",  matC[12], matC[13], matC[14], matC[15]);
       }
    printf("\n=== mncblas_sgemm ===\n\n");
    for (i = 0 ; i < NB_FOIS; i++) {
        start = _rdtsc () ;
        mncblas_sgemm(MNCblasRowMajor, MNCblasNoTrans,MNCblasNoTrans, VECSIZE, VECSIZE,VECSIZE, 1.0, mef,VECSIZE, mef2, VECSIZE,1.0, mef3, VECSIZE);
        end =_rdtsc() ;
        printf ("mncblas_sgemm %d : nombre de cycles: %Ld\n", i, end-start) ;
        calcul_flop ("sgemm ", VECSIZE * 3*VECSIZE *VECSIZE , end-start) ;
   }

   printf("\n=== mncblas_dgemm ===\n\n");
   for (i = 0 ; i < NB_FOIS; i++) {
       start = _rdtsc () ;
       mncblas_dgemm(MNCblasRowMajor, MNCblasNoTrans,MNCblasNoTrans, VECSIZE, VECSIZE,VECSIZE, 1.0, med,VECSIZE, med2, VECSIZE,1.0, med3, VECSIZE);
       end = _rdtsc () ;

       printf ("mncblas_dgemm %d : nombre de cycles: %Ld\n", i, end-start) ;
       calcul_flop ("dgemm ", VECSIZE * 3*VECSIZE *VECSIZE, end-start) ;
  }

  printf("\n=== mncblas_cgemm ===\n\n");
  complexe_float_t scal = (complexe_float_t){1.0,0.0};
  for (i = 0 ; i < NB_FOIS; i++) {
      start = _rdtsc () ;
        mncblas_cgemm(MNCblasRowMajor, MNCblasNoTrans,MNCblasNoTrans, VECSIZE, VECSIZE,VECSIZE, &scal, mecf,VECSIZE, mecf2, VECSIZE,&scal, mecf3, VECSIZE);
      end = _rdtsc () ;

      printf ("mncblas_cgemm %d : nombre de cycles: %Ld\n", i, end-start) ;
      calcul_flop ("cgemm ", VECSIZE * VECSIZE*23 * VECSIZE * 8, end-start) ;
 }

 printf("\n=== mncblas_zgemm ===\n\n");
 complexe_double_t sal = (complexe_double_t){1.0,0.0};
 for (i = 0 ; i < NB_FOIS; i++) {
     start = _rdtsc () ;
           mncblas_zgemm(MNCblasRowMajor, MNCblasNoTrans,MNCblasNoTrans, VECSIZE, VECSIZE,VECSIZE, &sal, mecd,VECSIZE, mecd2, VECSIZE,&sal, mecd3, VECSIZE);
     end = _rdtsc () ;

     printf ("mncblas_zgemm %d : nombre de cycles: %Ld\n", i, end-start) ;
     calcul_flop ("zgemm ", VECSIZE * VECSIZE*23 * VECSIZE * 8, end-start) ;
 }
}
