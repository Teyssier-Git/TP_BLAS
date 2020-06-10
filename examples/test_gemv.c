#include <stdio.h>
#include <x86intrin.h>

#include "mnblas.h"
#include "complexe.h"

#include "flop.h"

#define VECSIZE    6553

#define NB_FOIS    10

typedef float vfloat [VECSIZE];
typedef float mfloat [VECSIZE*VECSIZE];
typedef double vdouble [VECSIZE];
typedef double mdouble [VECSIZE*VECSIZE];
typedef complexe_float_t vCompFloat [VECSIZE];
typedef complexe_float_t mCompFloat [VECSIZE*VECSIZE];
typedef complexe_double_t vCompDouble [VECSIZE];
typedef complexe_double_t mCompDouble [VECSIZE*VECSIZE];

vfloat vef;
mfloat mef;
vdouble ved;
mdouble med;
vCompFloat vecf;
mCompFloat mecf;
vCompDouble vecd;
mCompDouble mecd;

void vector_init_s (vfloat V, float x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE; i++) {
        V [i] = x ;
    }
    return;
}

void matrix_init_s (mfloat M, float x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE*VECSIZE; i++) {
        M [i] = x ;
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

void matrix_init_d (mdouble M, double x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE*VECSIZE; i++) {
        M [i] = x ;
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

void matrix_init_c (mCompFloat M, complexe_float_t x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE*VECSIZE; i++) {
        M[i].real = x.real;
        M[i].imaginary = x.imaginary;
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

void matrix_init_z (mCompDouble M, complexe_double_t x) {
    register unsigned int i;
    for (i = 0; i < VECSIZE*VECSIZE; i++) {
        M[i].real = x.real;
        M[i].imaginary = x.imaginary;
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

void test_algorithme() {
    printf("TEST GEMV\n");

    printf("\n=== mncblas_sgemv ===\n\n");
    if (1) {
        float matA[4*4], vecX[4], vecY[4], vecYP[4];
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                matA[4*i+j] = i+1;
            }
            vecX[i] = i+1;
            vecY[i] = 2*(i+1);
            vecYP[i] = 2*(i+1);
        }
        mncblas_sgemv(MNCblasRowMajor, MNCblasNoTrans, 4, 4, 1.0, matA, 4, vecX, 1, 0.5, vecY, 1);
        printf("      | %.1f %.1f %.1f %.1f |   | %.1f |         | %.1f |   | %.1f |\n", vecX[0], matA[0], matA[1], matA[2], matA[3], vecYP[0], vecY[0]);
        printf("1.0 * | %.1f %.1f %.1f %.1f | * | %.1f | + 0.5 * | %.1f | = | %.1f |\n", vecX[1], matA[4], matA[5], matA[6], matA[7], vecYP[1], vecY[1]);
        printf("      | %.1f %.1f %.1f %.1f |   | %.1f |         | %.1f |   | %.1f |\n", vecX[2], matA[8], matA[9], matA[10], matA[11], vecYP[2], vecY[2]);
        printf("      | %.1f %.1f %.1f %.1f |   | %.1f |         | %.1f |   | %.1f |\n", vecX[3], matA[12], matA[13], matA[14], matA[15], vecYP[3], vecY[3]);
    }

    printf("\n=== mncblas_dgemv ===\n\n");
    if (1) {
        double matA[4*4], vecX[4], vecY[4], vecYP[4];
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                matA[4*i+j] = i+1;
            }
            vecX[i] = i+1;
            vecY[i] = 2*(i+1);
            vecYP[i] = 2*(i+1);
        }
        mncblas_dgemv(MNCblasRowMajor, MNCblasNoTrans, 4, 4, 1.0, matA, 4, vecX, 1, 0.5, vecY, 1);
        printf("      | %.1lf %.1lf %.1lf %.1lf |   | %.1lf |         | %.1lf |   | %.1lf |\n", vecX[0], matA[0], matA[1], matA[2], matA[3], vecYP[0], vecY[0]);
        printf("1.0 * | %.1lf %.1lf %.1lf %.1lf | * | %.1lf | + 0.5 * | %.1lf | = | %.1lf |\n", vecX[1], matA[4], matA[5], matA[6], matA[7], vecYP[1], vecY[1]);
        printf("      | %.1lf %.1lf %.1lf %.1lf |   | %.1lf |         | %.1lf |   | %.1lf |\n", vecX[2], matA[8], matA[9], matA[10], matA[11], vecYP[2], vecY[2]);
        printf("      | %.1lf %.1lf %.1lf %.1lf |   | %.1lf |         | %.1lf |   | %.1lf |\n", vecX[3], matA[12], matA[13], matA[14], matA[15], vecYP[3], vecY[3]);
    }

    printf("\n=== mncblas_cgemv ===\n\n");
    if (1) {
        complexe_float_t matA[4*4], vecX[4], rep[4], vecYP[4];
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                matA[4*i+j] = (complexe_float_t){i+1,0.0};
            }
            vecX[i] = (complexe_float_t){i+1,0.0};
            rep[i] = (complexe_float_t){2*(i+1),0.0};
            vecYP[i] = (complexe_float_t){2*(i+1),0.0};
        }
        complexe_float_t alpha = (complexe_float_t){1.0,0.0}, beta = (complexe_float_t){0.5,0.0};
        mncblas_cgemv(MNCblasRowMajor, MNCblasNoTrans, 4, 4, &alpha, matA, 4, vecX, 1, &beta, rep, 1);
        printf("      | %.1f %.1f %.1f %.1f |   | %.1f |         | %.1f |   | %.1f |\n", vecX[0].real, matA[0].real, matA[1].real, matA[2].real, matA[3].real, vecYP[0].real, rep[0].real);
        printf("1.0 * | %.1f %.1f %.1f %.1f | * | %.1f | + 0.5 * | %.1f | = | %.1f |\n", vecX[1].real, matA[4].real, matA[5].real, matA[6].real, matA[7].real, vecYP[1].real, rep[1].real);
        printf("      | %.1f %.1f %.1f %.1f |   | %.1f |         | %.1f |   | %.1f |\n", vecX[2].real, matA[8].real, matA[9].real, matA[10].real, matA[11].real, vecYP[2].real, rep[2].real);
        printf("      | %.1f %.1f %.1f %.1f |   | %.1f |         | %.1f |   | %.1f |\n", vecX[3].real, matA[12].real, matA[13].real, matA[14].real, matA[15].real, vecYP[3].real, rep[3].real);
    }

    printf("\n=== mncblas_zgemv ===\n\n");
    if (1) {
        complexe_double_t matA[4*4], vecX[4], rep[4], vecYP[4];
        for (int i=0; i<4; i++) {
            for (int j=0; j<4; j++) {
                matA[4*i+j] = (complexe_double_t){i+1,0.0};
            }
            vecX[i] = (complexe_double_t){i+1,0.0};
            rep[i] = (complexe_double_t){2*(i+1),0.0};
            vecYP[i] = (complexe_double_t){2*(i+1),0.0};
        }
        complexe_double_t alpha = (complexe_double_t){1.0,0.0}, beta = (complexe_double_t){0.5,0.0};
        mncblas_zgemv(MNCblasRowMajor, MNCblasNoTrans, 4, 4, &alpha, matA, 4, vecX, 1, &beta, rep, 1);
        printf("      | %.1lf %.1lf %.1lf %.1lf |   | %.1lf |         | %.1lf |   | %.1lf |\n", vecX[0].real, matA[0].real, matA[1].real, matA[2].real, matA[3].real, vecYP[0].real, rep[0].real);
        printf("1.0 * | %.1lf %.1lf %.1lf %.1lf | * | %.1lf | + 0.5 * | %.1lf | = | %.1lf |\n", vecX[1].real, matA[4].real, matA[5].real, matA[6].real, matA[7].real, vecYP[1].real, rep[1].real);
        printf("      | %.1lf %.1lf %.1lf %.1lf |   | %.1lf |         | %.1lf |   | %.1lf |\n", vecX[2].real, matA[8].real, matA[9].real, matA[10].real, matA[11].real, vecYP[2].real, rep[2].real);
        printf("      | %.1lf %.1lf %.1lf %.1lf |   | %.1lf |         | %.1lf |   | %.1lf |\n", vecX[3].real, matA[12].real, matA[13].real, matA[14].real, matA[15].real, vecYP[3].real, rep[3].real);
    }

}

void test_performance() {
    unsigned long long start, end;
    int i;

    init_flop () ;

    printf("TEST GEMV\n");

    printf("\n=== mncblas_sgemv ===\n\n");
    for (i = 0 ; i < NB_FOIS; i++) {
        start = _rdtsc () ;
            mncblas_sgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, VECSIZE, 1.0, mef, VECSIZE, vef, 1, 1.0, vef, 1);
        end =_rdtsc() ;
        printf ("mncblas_sgemv %d : nombre de cycles: %Ld\n", i, end-start) ;
        calcul_flop ("sgemv ", VECSIZE * (3 + VECSIZE * 2), end-start) ;
    }

    printf("\n=== mncblas_dgemv ===\n\n");
    for (i = 0 ; i < NB_FOIS; i++) {
        start = _rdtsc () ;
            mncblas_dgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, VECSIZE, 1.0, med, VECSIZE, ved, 1, 1.0, ved, 1) ;
        end = _rdtsc () ;

        printf ("mncblas_dgemv %d : nombre de cycles: %Ld\n", i, end-start) ;
        calcul_flop ("dgemv ", VECSIZE * (3 + VECSIZE * 2), end-start) ;
    }

    printf("\n=== mncblas_cgemv ===\n\n");
    complexe_float_t scal = (complexe_float_t){1.0,0.0};
    for (i = 0 ; i < NB_FOIS; i++) {
        start = _rdtsc () ;
            mncblas_cgemv(MNCblasRowMajor, MNCblasNoTrans, VECSIZE, VECSIZE, &scal, mecf, VECSIZE, vecf, 1, &scal, vecf, 1);
        end = _rdtsc () ;

        printf ("mncblas_cgemv %d : nombre de cycles: %Ld\n", i, end-start) ;
        calcul_flop ("cgemv ", VECSIZE * (14 + VECSIZE * 8), end-start) ;
    }

    printf("\n=== mncblas_zgemv ===\n\n");
    complexe_double_t sal = (complexe_double_t){1.0,0.0};
    for (i = 0 ; i < NB_FOIS; i++) {
        start = _rdtsc () ;
            mncblas_zgemv (MNCblasRowMajor, MNCblasNoTrans, VECSIZE, VECSIZE, &sal, mecd, VECSIZE, vecd, 1, &sal, vecd, 1) ;
        end = _rdtsc () ;

        printf ("mncblas_zgemv %d : nombre de cycles: %Ld\n", i, end-start) ;
        calcul_flop ("zgemv ", VECSIZE * (14 + VECSIZE * 8), end-start) ;
    }
}

int main (int argc, char **argv) {
    //test_algorithme();
    test_performance();
}
