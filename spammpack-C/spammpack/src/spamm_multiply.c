/** @file */

#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_SSE
#include <xmmintrin.h>
#endif

/* Some commonly used bit-patterns are:
 *
 * For 3D indices:
 *
 * 010010010010010010010010010010 = 0x12492492
 * 101101101101101101101101101101 = 0x2DB6DB6D
 *
 * For 2D indices:
 * 01010101010101010101010101010101 = 0x55555555
 * 10101010101010101010101010101010 = 0xaaaaaaaa
 */

/** 010010010010010010010010010010 = 0x12492492 */
#define MASK_3D_K  0x12492492

/** 101101101101101101101101101101 = 0x2DB6DB6D */
#define MASK_3D_IJ 0x2DB6DB6D

/** 01010101010101010101010101010101 = 0x55555555 */
#define MASK_2D_J  0x55555555

/** 10101010101010101010101010101010 = 0xaaaaaaaa */
#define MASK_2D_I  0xaaaaaaaa

/** @private k lookup table to speed up loop over indices. */
struct spamm_multiply_k_lookup_t
{
  /** The number of elements in the index array. */
  unsigned int size;

  /** An array of index values pointing into the linear index array for
   * different k values. */
  unsigned int *index;
};

/** @private Multiply a matrix by a scalar.
 *
 * @param alpha The scalar \f$\alpha\f$ that multiplies the matrix.
 * @param A The matrix.
 */
void
spamm_recursive_multiply_scalar (const float alpha,
    struct spamm_recursive_node_t *A,
    const unsigned int number_dimensions,
    const unsigned int tier,
    const unsigned int chunk_tier,
    const short use_linear_tree)
{
  unsigned int i;

  if (A == NULL) { return; }

  if (tier == chunk_tier)
  {
    A->norm2 = spamm_chunk_multiply_scalar(alpha, A->tree.chunk);
    A->norm = sqrt(A->norm2);
  }

  else
  {
    for (i = 0; i < ipow(2, number_dimensions); i++)
    {
      spamm_recursive_multiply_scalar(alpha, A->tree.child[i],
          number_dimensions, tier+1, chunk_tier, use_linear_tree);
    }
  }
}

/** Multiply two matrices, i.e. \f$ C = \alpha A \times B + \beta C\f$.
 *
 * @param tolerance The SpAMM tolerance of this product.
 * @param alpha The paramater \f$\alpha\f$.
 * @param A The matrix \f$A\f$.
 * @param B The matrix \f$B\f$.
 * @param beta The paramater \f$\beta\f$.
 * @param C The matrix \f$C\f$.
 * @param timer The timer to use.
 * @param kernel The stream kernel to use.
 *
 * @return The square of the Frobenius norm of the chunk.
 */
float
spamm_linear_multiply (const float tolerance,
    const float alpha,
    spamm_chunk_t *chunk_A,
    spamm_chunk_t *chunk_B,
    const float beta,
    spamm_chunk_t *chunk_C,
    struct spamm_timer_t *timer)
{
  unsigned int *index_A;
  unsigned int *index_B;

  unsigned int *stream;

  unsigned int N_contiguous;
  unsigned int index_length;

  unsigned int i, j;
  unsigned int stream_index;

  float *norm_A;
  float *norm_B;

  N_contiguous = spamm_chunk_get_N_contiguous(chunk_A);
  index_length = N_contiguous/SPAMM_N_KERNEL;

  index_A = malloc(ipow(index_length, 2)*sizeof(unsigned int));
  index_B = malloc(ipow(index_length, 2)*sizeof(unsigned int));

  for (i = 0; i < ipow(index_length, 2); i++)
  {
    index_A[i] = i;
    index_B[i] = i;
  }

  /* Sort indices along k index. */
  spamm_sort_masked(ipow(index_length, 2), index_A, MASK_2D_J);
  spamm_sort_masked(ipow(index_length, 2), index_B, MASK_2D_I);

  /* Sort within each k-block by descending norm. */
  norm_A = spamm_chunk_get_tier_norm(*spamm_chunk_get_number_tiers(chunk_A), chunk_A);
  norm_B = spamm_chunk_get_tier_norm(*spamm_chunk_get_number_tiers(chunk_B), chunk_B);

  for (i = 0; i < index_length; i++)
  {
    spamm_sort_norm(index_length, &index_A[i*index_length], norm_A);
    spamm_sort_norm(index_length, &index_B[i*index_length], norm_B);
  }

  /* Convolute. */
  stream = malloc(ipow(index_length, 2)*3*sizeof(unsigned int));

  for (i = 0, stream_index = 0; i < index_length; i++) {
    for (j = 0; j < index_length; j++)
    {
      if (norm_A[index_A[i]]*norm_B[index_B[j]] > tolerance)
      {
        stream[3*stream_index+0] = index_A[i];
        stream[3*stream_index+1] = index_B[i];
        stream[3*stream_index+2] = (index_A[i] & MASK_2D_I) | (index_B[j] & MASK_2D_J);
        stream_index++;
      }
    }
  }
  printf("[multiply] Added %u block products to stream\n", stream_index);

  /* Run kernel. */
  spamm_stream_kernel(stream_index, alpha, tolerance, stream, chunk_A, chunk_B, chunk_C);

  /* Free memory. */
  free(stream);
  free(index_A);
  free(index_B);

  return 0;
}

/** Multiply two matrices, i.e. \f$ C = \alpha A \times B + \beta C\f$.
 *
 * @param tolerance The SpAMM tolerance of this product.
 * @param alpha The paramater \f$\alpha\f$.
 * @param A The matrix \f$A\f$.
 * @param B The matrix \f$B\f$.
 * @param beta The paramater \f$\beta\f$.
 * @param C The matrix \f$C\f$.
 * @param timer The timer to use.
 * @param sgemm The external sgemm function to use.
 * @param tier The tier.
 * @param chunk_tier The contiguous tier.
 * @param use_linear_tree If set to 1 then we will switch to a linear tree at
 * chunk_tier.
 * @param number_products [out] The number of block products.
 */
void
spamm_recursive_multiply (const float tolerance,
    const float alpha,
    struct spamm_recursive_node_t *node_A,
    struct spamm_recursive_node_t *node_B,
    struct spamm_recursive_node_t **node_C,
    struct spamm_timer_t *timer,
    sgemm_func sgemm,
    const enum spamm_kernel_t kernel,
    const unsigned int number_dimensions_A,
    const unsigned int number_dimensions_B,
    const unsigned int number_dimensions_C,
    const unsigned int *const N,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    const unsigned int tier,
    const unsigned int chunk_tier,
    const short use_linear_tree,
    unsigned int *number_products)
{
  float beta = 1.0;
  unsigned int *new_N_lower;
  unsigned int *new_N_upper;

  short i, j, k;

  if (node_A == NULL || node_B == NULL) { return; }

  /* We have to allocate a new C block a tier up. */
  if (*node_C == NULL)
  {
    *node_C = spamm_recursive_new_node();
  }

  if (tier == chunk_tier)
  {
    (*node_C)->norm2 = spamm_chunk_multiply(tolerance, alpha,
        node_A->tree.chunk, node_B->tree.chunk, (*node_C)->tree.chunk, tier,
        chunk_tier, 0, 0, 0, timer, sgemm, kernel);
    (*node_C)->norm = sqrt((*node_C)->norm2);
  }

  else
  {
    if ((*node_C)->tree.child == NULL)
    {
      (*node_C)->tree.child = calloc(ipow(2, number_dimensions_C), sizeof(struct spamm_recursive_node_t*));
    }

    new_N_lower = calloc(number_dimensions_C, sizeof(unsigned int));
    new_N_upper = calloc(number_dimensions_C, sizeof(unsigned int));

    if (number_dimensions_A == 2 &&
        number_dimensions_B == 2 &&
        number_dimensions_C == 2)
    {
      for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
          for (k = 0; k < 2; k++)
          {
            new_N_lower[0] = N_lower[0]+(N_upper[0]-N_lower[0])/2*i;
            new_N_upper[0] = N_lower[0]+(N_upper[0]-N_lower[0])/2*(i+1);
            new_N_lower[1] = N_lower[1]+(N_upper[1]-N_lower[1])/2*j;
            new_N_upper[1] = N_lower[1]+(N_upper[1]-N_lower[1])/2*(j+1);

            if (node_A->norm*node_B->norm > tolerance)
            {
              spamm_recursive_multiply(tolerance, alpha,
                  node_A->tree.child[i+2*k], node_B->tree.child[k+2*j],
                  &(*node_C)->tree.child[i+2*j], timer, sgemm, kernel,
                  number_dimensions_A, number_dimensions_B,
                  number_dimensions_C, N, new_N_lower, new_N_upper, tier+1,
                  chunk_tier, use_linear_tree, number_products);
            }
          }
        }
      }
    }

    else
    {
      SPAMM_FATAL("not implemented\n");
    }

    free(new_N_lower);
    free(new_N_upper);
  }
}

/** Multiply two matrices, i.e. \f$ C = \alpha A \times B + \beta C\f$.
 *
 * @param tolerance The SpAMM tolerance of this product.
 * @param alpha The paramater \f$\alpha\f$.
 * @param A The matrix \f$A\f$.
 * @param B The matrix \f$B\f$.
 * @param beta The paramater \f$\beta\f$.
 * @param C The matrix \f$C\f$.
 * @param timer The timer to use.
 * @param sgemm The external sgemm function to use.
 */
void
spamm_multiply (const float tolerance,
    const float alpha,
    struct spamm_matrix_t *A,
    struct spamm_matrix_t *B,
    const float beta,
    struct spamm_matrix_t *C,
    struct spamm_timer_t *timer,
    sgemm_func sgemm,
    const enum spamm_kernel_t kernel,
    unsigned int *number_products)
{
  int dim;
  unsigned int *N_lower;
  unsigned int *N_upper;

  if (A->chunk_tier == 0)
  {
    spamm_chunk_multiply_scalar(beta, C->tree.chunk);
    spamm_chunk_multiply(tolerance, alpha, A->tree.chunk, B->tree.chunk,
        C->tree.chunk, 0, C->chunk_tier, 0, 0, 0, timer, sgemm, kernel);
  }

  else
  {
    spamm_recursive_multiply_scalar(beta, C->tree.recursive_tree,
        C->number_dimensions, 0, C->chunk_tier, C->use_linear_tree);

    N_lower = calloc(C->number_dimensions, sizeof(unsigned int));
    N_upper = calloc(C->number_dimensions, sizeof(unsigned int));

    for (dim = 0; dim < A->number_dimensions; dim++)
    {
      N_upper[dim] = A->N_padded;
    }

    spamm_recursive_multiply(tolerance, alpha, A->tree.recursive_tree,
        B->tree.recursive_tree, &(C->tree.recursive_tree), timer, sgemm,
        kernel, A->number_dimensions, B->number_dimensions,
        C->number_dimensions, A->N, N_lower, N_upper, 0, A->chunk_tier,
        A->use_linear_tree, number_products);

    free(N_lower);
    free(N_upper);
  }
}
