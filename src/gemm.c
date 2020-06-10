#include "../include/mnblas.h"
#include "../include/complexe.h"
#include <stdio.h>

void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
  MNCBLAS_TRANSPOSE TransB, const int M, const int N,
  const int K, const float alpha, const float *A,
  const int lda, const float *B, const int ldb,
  const float beta, float *C, const int ldc){
    #pragma omp parallel for
    for (int k = 0; k<M ; k++)
    {
      int K_M = k*M;
      #pragma omp parallel for
      for (int i = 0 ; i<M ; i++)
      {
        C [i + K_M] *= beta/alpha;
        #pragma omp parallel for
        for (int j = 0; j<N ; j++)
        {
          C [i + K_M] += A [j + K_M] * B [i + j*M];
        }
        C [i + K_M] *= alpha;
      }
    }

  }


  void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
    MNCBLAS_TRANSPOSE TransB, const int M, const int N,
    const int K, const double alpha, const double *A,
    const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc){
      #pragma omp parallel for
      for (int k = 0; k<M ; k++)
      {
        int K_M = k*M;
        #pragma omp parallel for
        for (int i = 0 ; i<M ; i++)
        {
          C [i + K_M] *= beta/alpha;
          #pragma omp parallel for
          for (int j = 0; j<N ; j++)
          {
            C [i + K_M] += A [j + K_M] * B [i + j*M];
          }
          C [i + K_M] *= alpha;
        }
      }

    }



    void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
      MNCBLAS_TRANSPOSE TransB, const int M, const int N,
      const int K, const void *alpha, const void *A,
      const int lda, const void *B, const int ldb,
      const void *beta, void *C, const int ldc){
        complexe_float_t *fA = (complexe_float_t *)A;
        complexe_float_t *fB = (complexe_float_t *)B;
        complexe_float_t *fC = (complexe_float_t *)C;
        complexe_float_t *al = (complexe_float_t *)alpha, *be = (complexe_float_t *)beta;
        complexe_float_t alp = *al, bet = *be;
        #pragma omp parallel for
        for (int k = 0; k<M ; k++)
        {
          int K_M = k*M;
          #pragma omp parallel for
          for (int i = 0 ; i<M ; i++)
          {
            fC [i + K_M] =mult_complexe_float(fC[i + K_M], div_complexe_float(bet, alp));
            #pragma omp parallel for
            for (int j = 0; j<N ; j++)
            {
              fC [i + K_M] =add_complexe_float(fC[i + K_M], mult_complexe_float(fA[j + K_M], fB[i + j*M]));
            }
            fC [i + K_M] = mult_complexe_float(fC[i + K_M], alp);
          }
        }

      }

      void mncblas_zgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
        MNCBLAS_TRANSPOSE TransB, const int M, const int N,
        const int K, const void *alpha, const void *A,
        const int lda, const void *B, const int ldb,
        const void *beta, void *C, const int ldc){
          complexe_double_t *dA = (complexe_double_t *)A;
          complexe_double_t *dB = (complexe_double_t *)B;
          complexe_double_t *dC = (complexe_double_t *)C;
          complexe_double_t *al = (complexe_double_t *)alpha, *be = (complexe_double_t *)beta;
          complexe_double_t alp = *al, bet = *be;
          #pragma omp parallel for
          for (int k = 0; k<M ; k++)
          {
            int K_M = k*M;
            #pragma omp parallel for
            for (int i = 0 ; i<M ; i++)
            {
              dC [i + K_M] =mult_complexe_double(dC[i + K_M], div_complexe_double(bet, alp));
              #pragma omp parallel for
              for (int j = 0; j<N ; j++)
              {
                dC [i + K_M] =add_complexe_double(dC[i + K_M], mult_complexe_double(dA[j + K_M], dB[i + j*M]));
              }
              dC [i + K_M] = mult_complexe_double(dC[i + K_M], alp);
            }
          }
        }
#include "../include/mnblas.h"
#include "../include/complexe.h"
#include <stdio.h>

void mncblas_sgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
  MNCBLAS_TRANSPOSE TransB, const int M, const int N,
  const int K, const float alpha, const float *A,
  const int lda, const float *B, const int ldb,
  const float beta, float *C, const int ldc){

    for (int k = 0; k<M ; k++)
                       {
                         int K_M = k*M;
                         for (int i = 0 ; i<M ; i++)
                            {
                             C [i + K_M] *= beta/alpha;
                             for (int j = 0; j<N ; j++)
                                 {
                                    C [i + K_M] += A [j + K_M] * B [i + j*M];
                                 }
                             C [i + K_M] *= alpha;
                           }
                        }

  }


  void mncblas_dgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
    MNCBLAS_TRANSPOSE TransB, const int M, const int N,
    const int K, const double alpha, const double *A,
    const int lda, const double *B, const int ldb,
    const double beta, double *C, const int ldc){
      for (int k = 0; k<M ; k++)
                         {
                           int K_M = k*M;
                           for (int i = 0 ; i<M ; i++)
                              {
                               C [i + K_M] *= beta/alpha;
                               for (int j = 0; j<N ; j++)
                                   {
                                      C [i + K_M] += A [j + K_M] * B [i + j*M];
                                   }
                               C [i + K_M] *= alpha;
                             }
                          }

    }



    void mncblas_cgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
      MNCBLAS_TRANSPOSE TransB, const int M, const int N,
      const int K, const void *alpha, const void *A,
      const int lda, const void *B, const int ldb,
      const void *beta, void *C, const int ldc){
        complexe_float_t *fA = (complexe_float_t *)A;
        complexe_float_t *fB = (complexe_float_t *)B;
        complexe_float_t *fC = (complexe_float_t *)C;
        complexe_float_t *al = (complexe_float_t *)alpha, *be = (complexe_float_t *)beta;
        complexe_float_t alp = *al, bet = *be;

        for (int k = 0; k<M ; k++)
                           {
                             int K_M = k*M;
                             for (int i = 0 ; i<M ; i++)
                                {
                                 fC [i + K_M] =mult_complexe_float(fC[i + K_M], div_complexe_float(bet, alp));
                                 for (int j = 0; j<N ; j++)
                                     {
                                        fC [i + K_M] =add_complexe_float(fC[i + K_M], mult_complexe_float(fA[j + K_M], fB[i + j*M]));
                                     }
                                 fC [i + K_M] = mult_complexe_float(fC[i + K_M], alp);
                               }
                            }

      }

      void mncblas_zgemm(MNCBLAS_LAYOUT layout, MNCBLAS_TRANSPOSE TransA,
        MNCBLAS_TRANSPOSE TransB, const int M, const int N,
        const int K, const void *alpha, const void *A,
        const int lda, const void *B, const int ldb,
        const void *beta, void *C, const int ldc){
          complexe_double_t *dA = (complexe_double_t *)A;
          complexe_double_t *dB = (complexe_double_t *)B;
          complexe_double_t *dC = (complexe_double_t *)C;
          complexe_double_t *al = (complexe_double_t *)alpha, *be = (complexe_double_t *)beta;
          complexe_double_t alp = *al, bet = *be;

          for (int k = 0; k<M ; k++)
                             {
                               int K_M = k*M;
                               for (int i = 0 ; i<M ; i++)
                                  {
                                   dC [i + K_M] =mult_complexe_double(dC[i + K_M], div_complexe_double(bet, alp));
                                   for (int j = 0; j<N ; j++)
                                       {
                                          dC [i + K_M] =add_complexe_double(dC[i + K_M], mult_complexe_double(dA[j + K_M], dB[i + j*M]));
                                       }
                                   dC [i + K_M] = mult_complexe_double(dC[i + K_M], alp);
                                 }
                              }
        }
