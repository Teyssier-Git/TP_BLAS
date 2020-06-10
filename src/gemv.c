#include "mnblas.h"
#include "complexe.h"
#include <stdio.h>

void mncblas_sgemv(const MNCBLAS_LAYOUT layout,
                 const MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *X, const int incX, const float beta,
                 float *Y, const int incY) {

    register unsigned int i;

    if (TransA == MNCblasNoTrans) {
        if (layout == MNCblasRowMajor) {
            #pragma omp parallel for
            for (i = 0; i < (M/incX+1) ; i ++) {
                float resLine = 0.0;
                for (int k=0; k<N; k++) {
                    resLine += A[N*i*incX+k]*X[k];
                }
                Y[i*incY] = alpha * resLine + beta * Y[i*incY];;
            }
        } else {
            #pragma omp parallel for
            for (i = 0; i < M/incX+1; i++) {
                float resLine = 0.0;
                for (int k=0; k<N; k++) {
                    resLine += A[M*k+i*incX]*X[k];
                }
                Y[i*incY] = alpha * resLine + beta * Y[i*incY];
            }
        }
    } else if (TransA == MNCblasTrans || TransA == MNCblasConjTrans) {
        if (layout == MNCblasRowMajor) {
            #pragma omp parallel for
            for (i = 0; i < N/incX+1; i++) {
                float resLine = 0.0;
                for (int k=0; k<M; k++) {
                    resLine += A[N*k+i*incX]*X[k];
                }
                Y[i*incY] = alpha * resLine + beta * Y[i*incY];
            }
        } else {
            #pragma omp parallel for
            for (i = 0; i < N/incX+1 ; i++) {
                float resLine = 0.0;
                for (int k=0; k<M; k++) {
                    resLine += A[M*i*incX+k]*X[k];
                }
                Y[i*incY] = alpha * resLine + beta * Y[i*incY];
            }
        }
    }
}


void mncblas_dgemv(MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY) {

    register unsigned int i;

    if (TransA == MNCblasNoTrans) {
        if (layout == MNCblasRowMajor) {
            #pragma omp parallel for
            for (i = 0; i < M/incX+1 ; i++) {
                double resLine = 0.0;
                for (int k=0; k<N; k++) {
                    resLine += A[N*i*incX+k]*X[k];
                }
                Y[i*incY] = alpha * resLine + beta * Y[i*incY];
            }
        } else {
            #pragma omp parallel for
            for (i = 0; i < M /incX+1; i++) {
                double resLine = 0.0;
                for (int k=0; k<N; k++) {
                    resLine += A[M*k+i*incX]*X[k];
                }
                Y[i*incY] = alpha * resLine + beta * Y[i*incY];
            }
        }
    } else if (TransA == MNCblasTrans || TransA == MNCblasConjTrans) {
        if (layout == MNCblasRowMajor) {
            #pragma omp parallel for
            for (i = 0; i < N/incX+1 ; i++) {
                double resLine = 0.0;
                for (int k=0; k<M; k++) {
                    resLine += A[N*k+i*incX]*X[k];
                }
                Y[i*incY] = alpha * resLine + beta * Y[i*incY];
            }
        } else {
            #pragma omp parallel for
            for (i = 0; i < N/incX+1 ; i++) {
                double resLine = 0.0;
                for (int k=0; k<M; k++) {
                    resLine += A[M*i*incX+k]*X[k];
                }
                Y[i*incY] = alpha * resLine + beta * Y[i*incY];
            }
        }
    }
}


void mncblas_cgemv(MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY) {

    register unsigned int i;

    complexe_float_t *x = (complexe_float_t *)X, *y = (complexe_float_t *)Y, *a = (complexe_float_t *)A;
    complexe_float_t *al = (complexe_float_t *)alpha, *be = (complexe_float_t *)beta;
    complexe_float_t alp = *al, bet = *be;

    if (TransA == MNCblasNoTrans) {
        if (layout == MNCblasRowMajor) {
            #pragma omp parallel for
            for (i = 0; i < M/incX+1 ; i++) {
                complexe_float_t resLine = (complexe_float_t){0.0,0.0};
                for (int k=0; k<N; k++) {
                    resLine = add_complexe_float(resLine,mult_complexe_float(a[N*i*incX+k],x[k]));
                }
                y[i*incY] = add_complexe_float(mult_complexe_float(alp, resLine), mult_complexe_float(bet, y[i*incY]));
            }
        } else {
            #pragma omp parallel for
            for (i = 0; i < M/incX+1 ; i++) {
                complexe_float_t resLine = (complexe_float_t){0.0,0.0};
                for (int k=0; k<N; k++) {
                    resLine = add_complexe_float(resLine, mult_complexe_float(a[M*k+i*incX], x[k]));
                }
                y[i*incY] = add_complexe_float(mult_complexe_float(alp, resLine), mult_complexe_float(bet, y[i*incY]));
            }
        }
    } else if (TransA == MNCblasTrans) {
        if (layout == MNCblasRowMajor) {
            #pragma omp parallel for
            for (i = 0; i < N/incX+1 ; i++) {
                complexe_float_t resLine = (complexe_float_t){0.0,0.0};
                for (int k=0; k<M; k++) {
                    resLine = add_complexe_float(resLine, mult_complexe_float(a[N*k+i*incX], x[k]));
                }
                y[i*incY] = add_complexe_float(mult_complexe_float(alp, resLine), mult_complexe_float(bet, y[i*incY]));
            }
        } else {
            #pragma omp parallel for
            for (i = 0; i < N; i++) {
                complexe_float_t resLine = (complexe_float_t){0.0,0.0};
                for (int k=0; k<M; k++) {
                    resLine = add_complexe_float(resLine, mult_complexe_float(a[M*i*incX+k], x[k]));
                }
                y[i*incY] = add_complexe_float(mult_complexe_float(alp, resLine), mult_complexe_float(bet, y[i*incY]));
            }
        }
    } else {
        if (layout == MNCblasRowMajor) {
            #pragma omp parallel for
            for (i = 0; i < N/incX+1 ; i++) {
                complexe_float_t resLine = (complexe_float_t){0.0,0.0};
                for (int k=0; k<M; k++) {
                    resLine = add_complexe_float(resLine, mult_complexe_float(conj_complexe_float(a[N*k+i*incX]), x[k]));
                }
                y[i*incY] = add_complexe_float(mult_complexe_float(alp, resLine), mult_complexe_float(bet, y[i*incY]));
            }
        } else {
            #pragma omp parallel for
            for (i = 0; i < N/incX+1 ; i++) {
                complexe_float_t resLine = (complexe_float_t){0.0,0.0};
                for (int k=0; k<M; k++) {
                    resLine = add_complexe_float(resLine, mult_complexe_float(conj_complexe_float(a[M*i*incX+k]), x[k]));
                }
                y[i*incY] = add_complexe_float(mult_complexe_float(alp, resLine), mult_complexe_float(bet, y[i*incY]));
            }
        }
    }
}

void mncblas_zgemv(MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY) {

    register unsigned int i;

    complexe_double_t *x = (complexe_double_t *)X, *y = (complexe_double_t *)Y, *a = (complexe_double_t *)A;
    complexe_double_t *al = (complexe_double_t *)alpha, *be = (complexe_double_t *)beta;
    complexe_double_t alp = *al, bet = *be;

    if (TransA == MNCblasNoTrans) {
        if (layout == MNCblasRowMajor) {
            #pragma omp parallel for
            for (i = 0; i < M/incX+1 ; i++) {
                complexe_double_t resLine = (complexe_double_t){0.0,0.0};
                for (int k=0; k<N; k++) {
                    resLine = add_complexe_double(resLine,mult_complexe_double(a[N*i*incX+k],x[k]));
                }
                y[i*incY] = add_complexe_double(mult_complexe_double(alp, resLine), mult_complexe_double(bet, y[i*incY]));
            }
        } else {
            #pragma omp parallel for
            for (i = 0; i < M/incX+1 ; i++) {
                complexe_double_t resLine = (complexe_double_t){0.0,0.0};
                for (int k=0; k<N; k++) {
                    resLine = add_complexe_double(resLine, mult_complexe_double(a[M*k+i*incX], x[k]));
                }
                y[i*incY] = add_complexe_double(mult_complexe_double(alp, resLine), mult_complexe_double(bet, y[i*incY]));
            }
        }
    } else if (TransA == MNCblasTrans) {
        if (layout == MNCblasRowMajor) {
            #pragma omp parallel for
            for (i = 0; i < N/incX+1 ; i++) {
                complexe_double_t resLine = (complexe_double_t){0.0,0.0};
                for (int k=0; k<M; k++) {
                    resLine = add_complexe_double(resLine, mult_complexe_double(a[N*k+i*incX], x[k]));
                }
                y[i*incY] = add_complexe_double(mult_complexe_double(alp, resLine), mult_complexe_double(bet, y[i*incY]));
            }
        } else {
            #pragma omp parallel for
            for (i = 0; i < N ; i++) {
                complexe_double_t resLine = (complexe_double_t){0.0,0.0};
                for (int k=0; k<M; k++) {
                    resLine = add_complexe_double(resLine, mult_complexe_double(a[M*i*incX+k], x[k]));
                }
                y[i*incY] = add_complexe_double(mult_complexe_double(alp, resLine), mult_complexe_double(bet, y[i*incY]));
            }
        }
    } else {
        if (layout == MNCblasRowMajor) {
            #pragma omp parallel for
            for (i = 0; i < N/incX+1 ; i++) {
                complexe_double_t resLine = (complexe_double_t){0.0,0.0};
                for (int k=0; k<M; k++) {
                    resLine = add_complexe_double(resLine, mult_complexe_double(conj_complexe_double(a[N*k+i*incX]), x[k]));
                }
                y[incY*i] = add_complexe_double(mult_complexe_double(alp, resLine), mult_complexe_double(bet, y[i*incY]));
            }
        } else {
            #pragma omp parallel for
            for (i = 0; i < N ; i++) {
                complexe_double_t resLine = (complexe_double_t){0.0,0.0};
                for (int k=0; k<M; k++) {
                    resLine = add_complexe_double(resLine, mult_complexe_double(conj_complexe_double(a[M*i*incX+k]), x[k]));
                }
                y[i*incY] = add_complexe_double(mult_complexe_double(alp, resLine), mult_complexe_double(bet, y[i*incY]));
            }
        }
    }
}
