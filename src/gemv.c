#include "mnblas.h"
#include "complexe.h"
#include <stdio.h>

void mncblas_sgemv(const MNCBLAS_LAYOUT layout,
                 const MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const float alpha, const float *A, const int lda,
                 const float *X, const int incX, const float beta,
                 float *Y, const int incY) {

    register unsigned int i, j;

    if (TransA == MNCblasNoTrans) {
        if (layout == MNCblasRowMajor) {
            for (i = 0, j = 0 ; i < M ; i += incX, j += incY) {
                float resLine = 0.0;
                for (int k=0; k<N; k++) {
                    resLine += A[N*i+k]*X[k];
                }
                Y[j] = alpha * resLine + beta * Y[j];
            }
        } else {
            for (i = 0, j = 0 ; i < M ; i += incX, j += incY) {
                float resLine = 0.0;
                for (int k=0; k<N; k++) {
                    resLine += A[M*k+i]*X[k];
                }
                Y[j] = alpha * resLine + beta * Y[j];
            }
        }
    } else if (TransA == MNCblasTrans || TransA == MNCblasConjTrans) {
        if (layout == MNCblasRowMajor) {
            for (i = 0, j = 0 ; i < N ; i += incX, j += incY) {
                float resLine = 0.0;
                for (int k=0; k<M; k++) {
                    resLine += A[N*k+i]*X[k];
                }
                Y[j] = alpha * resLine + beta * Y[j];
            }
        } else {
            for (i = 0, j = 0 ; i < N ; i += incX, j += incY) {
                float resLine = 0.0;
                for (int k=0; k<M; k++) {
                    resLine += A[M*i+k]*X[k];
                }
                Y[j] = alpha * resLine + beta * Y[j];
            }
        }
    }
}


void mncblas_dgemv(MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX, const double beta,
                 double *Y, const int incY) {

    register unsigned int i, j;

    if (TransA == MNCblasNoTrans) {
        if (layout == MNCblasRowMajor) {
            for (i = 0, j = 0 ; i < M ; i += incX, j += incY) {
                double resLine = 0.0;
                for (int k=0; k<N; k++) {
                    resLine += A[N*i+k]*X[k];
                }
                Y[j] = alpha * resLine + beta * Y[j];
            }
        } else {
            for (i = 0, j = 0 ; i < M ; i += incX, j += incY) {
                double resLine = 0.0;
                for (int k=0; k<N; k++) {
                    resLine += A[M*k+i]*X[k];
                }
                Y[j] = alpha * resLine + beta * Y[j];
            }
        }
    } else if (TransA == MNCblasTrans || TransA == MNCblasConjTrans) {
        if (layout == MNCblasRowMajor) {
            for (i = 0, j = 0 ; i < N ; i += incX, j += incY) {
                double resLine = 0.0;
                for (int k=0; k<M; k++) {
                    resLine += A[N*k+i]*X[k];
                }
                Y[j] = alpha * resLine + beta * Y[j];
            }
        } else {
            for (i = 0, j = 0 ; i < N ; i += incX, j += incY) {
                double resLine = 0.0;
                for (int k=0; k<M; k++) {
                    resLine += A[M*i+k]*X[k];
                }
                Y[j] = alpha * resLine + beta * Y[j];
            }
        }
    }
}


void mncblas_cgemv(MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY) {

    register unsigned int i, j;

    complexe_float_t *x = (complexe_float_t *)X, *y = (complexe_float_t *)Y, *a = (complexe_float_t *)A;
    complexe_float_t *al = (complexe_float_t *)alpha, *be = (complexe_float_t *)beta;
    complexe_float_t alp = *al, bet = *be;

    if (TransA == MNCblasNoTrans) {
        if (layout == MNCblasRowMajor) {
            for (i = 0, j = 0 ; i < M ; i += incX, j += incY) {
                complexe_float_t resLine = (complexe_float_t){0.0,0.0};
                for (int k=0; k<N; k++) {
                    resLine = add_complexe_float(resLine,mult_complexe_float(a[N*i+k],x[k]));
                }
                y[j] = add_complexe_float(mult_complexe_float(alp, resLine), mult_complexe_float(bet, y[j]));
            }
        } else {
            for (i = 0, j = 0 ; i < M ; i += incX, j += incY) {
                complexe_float_t resLine = (complexe_float_t){0.0,0.0};
                for (int k=0; k<N; k++) {
                    resLine = add_complexe_float(resLine, mult_complexe_float(a[M*k+i], x[k]));
                }
                y[j] = add_complexe_float(mult_complexe_float(alp, resLine), mult_complexe_float(bet, y[j]));
            }
        }
    } else if (TransA == MNCblasTrans) {
        if (layout == MNCblasRowMajor) {
            for (i = 0, j = 0 ; i < N ; i += incX, j += incY) {
                complexe_float_t resLine = (complexe_float_t){0.0,0.0};
                for (int k=0; k<M; k++) {
                    resLine = add_complexe_float(resLine, mult_complexe_float(a[N*k+i], x[k]));
                }
                y[j] = add_complexe_float(mult_complexe_float(alp, resLine), mult_complexe_float(bet, y[j]));
            }
        } else {
            for (i = 0, j = 0 ; i < N ; i += incX, j += incY) {
                complexe_float_t resLine = (complexe_float_t){0.0,0.0};
                for (int k=0; k<M; k++) {
                    resLine = add_complexe_float(resLine, mult_complexe_float(a[M*i+k], x[k]));
                }
                y[j] = add_complexe_float(mult_complexe_float(alp, resLine), mult_complexe_float(bet, y[j]));
            }
        }
    } else {
        if (layout == MNCblasRowMajor) {
            for (i = 0, j = 0 ; i < N ; i += incX, j += incY) {
                complexe_float_t resLine = (complexe_float_t){0.0,0.0};
                for (int k=0; k<M; k++) {
                    resLine = add_complexe_float(resLine, mult_complexe_float(conj_complexe_float(a[N*k+i]), x[k]));
                }
                y[j] = add_complexe_float(mult_complexe_float(alp, resLine), mult_complexe_float(bet, y[j]));
            }
        } else {
            for (i = 0, j = 0 ; i < N ; i += incX, j += incY) {
                complexe_float_t resLine = (complexe_float_t){0.0,0.0};
                for (int k=0; k<M; k++) {
                    resLine = add_complexe_float(resLine, mult_complexe_float(conj_complexe_float(a[M*i+k]), x[k]));
                }
                y[j] = add_complexe_float(mult_complexe_float(alp, resLine), mult_complexe_float(bet, y[j]));
            }
        }
    }
}

void mncblas_zgemv(MNCBLAS_LAYOUT layout,
                 MNCBLAS_TRANSPOSE TransA, const int M, const int N,
                 const void *alpha, const void *A, const int lda,
                 const void *X, const int incX, const void *beta,
                 void *Y, const int incY) {

    register unsigned int i, j;

    complexe_double_t *x = (complexe_double_t *)X, *y = (complexe_double_t *)Y, *a = (complexe_double_t *)A;
    complexe_double_t *al = (complexe_double_t *)alpha, *be = (complexe_double_t *)beta;
    complexe_double_t alp = *al, bet = *be;

    if (TransA == MNCblasNoTrans) {
        if (layout == MNCblasRowMajor) {
            for (i = 0, j = 0 ; i < M ; i += incX, j += incY) {
                complexe_double_t resLine = (complexe_double_t){0.0,0.0};
                for (int k=0; k<N; k++) {
                    resLine = add_complexe_double(resLine,mult_complexe_double(a[N*i+k],x[k]));
                }
                y[j] = add_complexe_double(mult_complexe_double(alp, resLine), mult_complexe_double(bet, y[j]));
            }
        } else {
            for (i = 0, j = 0 ; i < M ; i += incX, j += incY) {
                complexe_double_t resLine = (complexe_double_t){0.0,0.0};
                for (int k=0; k<N; k++) {
                    resLine = add_complexe_double(resLine, mult_complexe_double(a[M*k+i], x[k]));
                }
                y[j] = add_complexe_double(mult_complexe_double(alp, resLine), mult_complexe_double(bet, y[j]));
            }
        }
    } else if (TransA == MNCblasTrans) {
        if (layout == MNCblasRowMajor) {
            for (i = 0, j = 0 ; i < N ; i += incX, j += incY) {
                complexe_double_t resLine = (complexe_double_t){0.0,0.0};
                for (int k=0; k<M; k++) {
                    resLine = add_complexe_double(resLine, mult_complexe_double(a[N*k+i], x[k]));
                }
                y[j] = add_complexe_double(mult_complexe_double(alp, resLine), mult_complexe_double(bet, y[j]));
            }
        } else {
            for (i = 0, j = 0 ; i < N ; i += incX, j += incY) {
                complexe_double_t resLine = (complexe_double_t){0.0,0.0};
                for (int k=0; k<M; k++) {
                    resLine = add_complexe_double(resLine, mult_complexe_double(a[M*i+k], x[k]));
                }
                y[j] = add_complexe_double(mult_complexe_double(alp, resLine), mult_complexe_double(bet, y[j]));
            }
        }
    } else {
        if (layout == MNCblasRowMajor) {
            for (i = 0, j = 0 ; i < N ; i += incX, j += incY) {
                complexe_double_t resLine = (complexe_double_t){0.0,0.0};
                for (int k=0; k<M; k++) {
                    resLine = add_complexe_double(resLine, mult_complexe_double(conj_complexe_double(a[N*k+i]), x[k]));
                }
                y[j] = add_complexe_double(mult_complexe_double(alp, resLine), mult_complexe_double(bet, y[j]));
            }
        } else {
            for (i = 0, j = 0 ; i < N ; i += incX, j += incY) {
                complexe_double_t resLine = (complexe_double_t){0.0,0.0};
                for (int k=0; k<M; k++) {
                    resLine = add_complexe_double(resLine, mult_complexe_double(conj_complexe_double(a[M*i+k]), x[k]));
                }
                y[j] = add_complexe_double(mult_complexe_double(alp, resLine), mult_complexe_double(bet, y[j]));
            }
        }
    }
}
