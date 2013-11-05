#ifndef __SPAMM_BLAS_H
#define __SPAMM_BLAS_H

#ifdef __cplusplus
#define __BEGIN_DECLARATIONS extern "C" {
#define __END_DECLARATIONS }
#else
#define __BEGIN_DECLARATIONS
#define __END_DECLARATIONS
#endif

__BEGIN_DECLARATIONS

#ifdef ADD_SGEMM_EXTERNAL_DECLARATION
void sgemm_ (char *transA, char *transB, int *M, int *N, int *K, float *alpha,
    float *A, int *LDA, float *B, int *LDB, float *beta, float *C, int *LDC);
#endif

#ifdef ADD_DGEMM_EXTERNAL_DECLARATION
void dgemm_ (char *transA, char *transB, int *M, int *N, int *K, double *alpha,
    double *A, int *LDA, double *B, int *LDB, double *beta, double *C, int *LDC);
#endif

__END_DECLARATIONS

#endif
