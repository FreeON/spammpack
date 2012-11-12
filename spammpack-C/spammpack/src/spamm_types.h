/** @file */

#ifndef __SPAMM_TYPES_H
#define __SPAMM_TYPES_H

typedef void (*sgemm_func) (char *transA, char *transB, int *M, int *N, int
    *K, float *alpha, float *A, int *LDA, float *B, int *LDB, float *beta,
    float *C, int *LDC);

#ifdef ADD_SGEMM_EXTERNAL_DECLARATION
void sgemm_ (char *transA, char *transB, int *M, int *N, int *K, float *alpha,
    float *A, int *LDA, float *B, int *LDB, float *beta, float *C, int *LDC);
#endif

/** The layout type for the layout of the basic matrix blocks on the kernel
 * tier.
 */
enum spamm_layout_t
{
  /** Layout in row-major order. */
  row_major,

  /** Layout in column-major order. */
  column_major,

  /** Layout in Z-curve order. */
  Z_curve,

  /** Layout as a dense kernel block, in row-major order. */
  dense_column_major
};

struct spamm_multiply_stream_t;
struct spamm_matrix_t;
struct spamm_hashed_t;
struct spamm_hashed_node_t;
struct spamm_hashed_data_t;
struct spamm_recursive_t;
struct spamm_recursive_node_t;

/** A data chunk. */
typedef void spamm_chunk_t;

#endif
