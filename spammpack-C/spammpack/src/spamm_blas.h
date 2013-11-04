#ifndef __SPAMM_BLAS_H
#define __SPAMM_BLAS_H

#ifdef ADD_SGEMM_EXTERNAL_DECLARATION
void sgemm_ (char *transA, char *transB, int *M, int *N, int *K, float *alpha,
    float *A, int *LDA, float *B, int *LDB, float *beta, float *C, int *LDC);
#endif

#ifdef ADD_DGEMM_EXTERNAL_DECLARATION
void dgemm_ (char *transA, char *transB, int *M, int *N, int *K, double *alpha,
    double *A, int *LDA, double *B, int *LDB, double *beta, double *C, int *LDC);
#endif


#endif
