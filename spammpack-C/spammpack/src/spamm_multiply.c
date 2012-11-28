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

  if(A == NULL) { return; }

  if(tier == chunk_tier)
  {
    A->norm2 = spamm_chunk_multiply_scalar(alpha, A->tree.chunk);
    A->norm = sqrt(A->norm2);
  }

  else
  {
    for(i = 0; i < ipow(2, number_dimensions); i++)
    {
      spamm_recursive_multiply_scalar(alpha, A->tree.child[i],
          number_dimensions, tier+1, chunk_tier, use_linear_tree);
    }
  }
}

#define SPAMM_MULTIPLY_DEBUG
#define RUN_ASSEMBLY_KERNEL

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
    spamm_chunk_t *chunk_C)
{
  unsigned int *index_A;
  unsigned int *index_B;

  unsigned int *stream;

  unsigned int N_contiguous;
  unsigned int index_length;

  unsigned int i;
  unsigned int j_A, j_B;
  unsigned int stream_index;

  float *norm_A;
  float *norm_B;

  N_contiguous = spamm_chunk_get_N_contiguous(chunk_A);
  index_length = N_contiguous/SPAMM_N_KERNEL;

  index_A = malloc(ipow(index_length, 2)*sizeof(unsigned int));
  index_B = malloc(ipow(index_length, 2)*sizeof(unsigned int));

  for(i = 0; i < ipow(index_length, 2); i++)
  {
    index_A[i] = i;
    index_B[i] = i;
  }

  /* Sort indices along k index. */
  spamm_sort_masked(ipow(index_length, 2), index_A, MASK_2D_J);
  spamm_sort_masked(ipow(index_length, 2), index_B, MASK_2D_I);

  /* Sort within each k-block by descending norm. */
  norm_A = spamm_chunk_get_tier_norm(*spamm_chunk_get_number_tiers(chunk_A)-SPAMM_KERNEL_DEPTH-1, chunk_A);
  norm_B = spamm_chunk_get_tier_norm(*spamm_chunk_get_number_tiers(chunk_B)-SPAMM_KERNEL_DEPTH-1, chunk_B);

  for(i = 0; i < index_length; i++)
  {
    spamm_sort_norm(index_length, &index_A[i*index_length], norm_A);
    spamm_sort_norm(index_length, &index_B[i*index_length], norm_B);
  }

#ifdef SPAMM_MULTIPLY_DEBUG
  printf("potentially %u products\n", ipow(index_length, 3));
#endif

  /* Convolute. */
  stream = malloc(ipow(index_length, 3)*3*sizeof(unsigned int));

#ifdef SPAMM_MULTIPLY_DEBUG
  printf("stream (%p):\n", stream);
#endif
  for(i = 0, stream_index = 0; i < index_length; i++)
  {
    for(j_A = i*index_length; j_A < (i+1)*index_length; j_A++) {
      for(j_B = i*index_length; j_B < (i+1)*index_length; j_B++)
      {
#ifdef SPAMM_MULTIPLY_DEBUG
        printf("comparing norms: %e (norm_A[%u]) * %e (norm_B[%u]) = %e",
            norm_A[index_A[j_A]],
            index_A[j_A],
            norm_B[index_B[j_B]],
            index_B[j_B],
            norm_A[index_A[j_A]]*norm_B[index_B[j_B]]);
#endif
        if(norm_A[index_A[j_A]]*norm_B[index_B[j_B]] > tolerance)
        {
          stream[3*stream_index+0] = index_A[j_A];
          stream[3*stream_index+1] = index_B[j_B];
          stream[3*stream_index+2] = (index_A[j_A] & MASK_2D_I) | (index_B[j_B] & MASK_2D_J);
#ifdef SPAMM_MULTIPLY_DEBUG
          printf(" -> adding stream[%u] = { %u, %u, %u }\n", stream_index,
              stream[3*stream_index+0],
              stream[3*stream_index+1],
              stream[3*stream_index+2]);
#endif
          stream_index++;
        }

        else
        {
#ifdef SPAMM_MULTIPLY_DEBUG
          printf(", done\n");
#endif
          break;
        }
      }
    }
  }

#ifdef SPAMM_MULTIPLY_DEBUG
  printf("[multiply] Added %u (out of %u possible) block products to stream\n", stream_index, ipow(index_length, 3));
#endif

  /* Run kernel. */
#ifdef RUN_ASSEMBLY_KERNEL
  spamm_stream_kernel(stream_index, alpha, tolerance, stream, chunk_A, chunk_B, chunk_C);
#else
  norm_A = spamm_chunk_get_tier_norm(*spamm_chunk_get_number_tiers(chunk_A)-1, chunk_A);
  norm_B = spamm_chunk_get_tier_norm(*spamm_chunk_get_number_tiers(chunk_B)-1, chunk_B);

  norm_C = spamm_chunk_get_tier_norm(*spamm_chunk_get_number_tiers(chunk_C)-1, chunk_C);
  norm2_C = spamm_chunk_get_tier_norm2(*spamm_chunk_get_number_tiers(chunk_C)-1, chunk_C);

  matrix_A = spamm_chunk_get_matrix(chunk_A);
  matrix_B = spamm_chunk_get_matrix(chunk_B);
  matrix_C = spamm_chunk_get_matrix(chunk_C);

#ifdef SPAMM_MULTIPLY_DEBUG
  printf("A: "); spamm_print_chunk(chunk_A);
  printf("B: "); spamm_print_chunk(chunk_B);
  printf("C: "); spamm_print_chunk(chunk_C);
#endif

#ifdef SPAMM_MULTIPLY_DEBUG
  printf("starting to calculate product\n");
  fflush(stdout);
#endif
  for(i_stream = 0; i_stream < stream_index; i_stream++)
  {
#ifdef SPAMM_MULTIPLY_DEBUG
    printf("lin.index (%u,%u,%u)\n", stream[3*i_stream+0], stream[3*i_stream+1], stream[3*i_stream+2]);
#endif

    for(i = 0; i < SPAMM_N_KERNEL_BLOCKED; i++) {
      for(j = 0; j < SPAMM_N_KERNEL_BLOCKED; j++) {
        for(k = 0; k < SPAMM_N_KERNEL_BLOCKED; k++)
        {
          norm_offset_A = stream[3*i_stream+0]*SPAMM_N_KERNEL_BLOCKED*SPAMM_N_KERNEL_BLOCKED
            +i*SPAMM_N_KERNEL_BLOCKED+k;
          norm_offset_B = stream[3*i_stream+1]*SPAMM_N_KERNEL_BLOCKED*SPAMM_N_KERNEL_BLOCKED
            +k*SPAMM_N_KERNEL_BLOCKED+j;
          norm_offset_C = stream[3*i_stream+2]*SPAMM_N_KERNEL_BLOCKED*SPAMM_N_KERNEL_BLOCKED
            +i*SPAMM_N_KERNEL_BLOCKED+j;

          offset_A = stream[3*i_stream+0]*SPAMM_N_KERNEL*SPAMM_N_KERNEL
            +(i*SPAMM_N_KERNEL_BLOCKED+k)*SPAMM_N_BLOCK*SPAMM_N_BLOCK;
          offset_B = stream[3*i_stream+1]*SPAMM_N_KERNEL*SPAMM_N_KERNEL
            +(k*SPAMM_N_KERNEL_BLOCKED+j)*SPAMM_N_BLOCK*SPAMM_N_BLOCK;
          offset_C = stream[3*i_stream+2]*SPAMM_N_KERNEL*SPAMM_N_KERNEL
            +(i*SPAMM_N_KERNEL_BLOCKED+j)*SPAMM_N_BLOCK*SPAMM_N_BLOCK;

#ifdef SPAMM_MULTIPLY_DEBUG
          printf("(%u,%u,%u) -> norm_offset_C = %u, offset_C = %u\n", i, j, k, norm_offset_C, offset_C);
#endif

          if(norm_A[norm_offset_A]*norm_B[norm_offset_B] > tolerance)
          {
            for(i_block = 0; i_block < SPAMM_N_BLOCK; i_block++) {
              for(j_block = 0; j_block < SPAMM_N_BLOCK; j_block++) {
                for(k_block = 0; k_block < SPAMM_N_BLOCK; k_block++)
                {
                  matrix_C[offset_C+i_block*SPAMM_N_BLOCK+j_block] +=
                    alpha
                    *matrix_A[offset_A+i_block*SPAMM_N_BLOCK+k_block]
                    *matrix_B[offset_B+k_block*SPAMM_N_BLOCK+j_block];
                }

                /* Update C norm. */
#ifdef SPAMM_MULTIPLY_DEBUG
                printf("(%u,%u) -> row_maj = %u\n", i_block, j_block, i_block*SPAMM_N_BLOCK+j_block);
                fflush(stdout);
#endif
              }
            }
          }
        }
      }
    }
  }
#endif

  /* Update norms. */
  SPAMM_WARN("FIXME: update norms on product\n");

  /* Free memory. */
  free(stream);
  free(index_A);
  free(index_B);

  return 0;
}

/** Multiply two SpAMM chunks.
 */
float
spamm_chunk_multiply (const float tolerance,
    const float alpha,
    spamm_chunk_t *chunk_A,
    spamm_chunk_t *chunk_B,
    spamm_chunk_t *chunk_C,
    sgemm_func sgemm)
{
  unsigned int i, j, k;

  float *norm_A;
  float *norm_B;

  float *norm_C;
  float *norm2_C;

  float alpha_sgemm = alpha;
  float beta = 1.0;

  int N_contiguous;

  float *matrix_A;
  float *matrix_B;
  float *matrix_C;

  short use_linear_tree;

  unsigned int number_dimensions;

  if(chunk_A == NULL || chunk_B == NULL) { return 0.0; }

  use_linear_tree = *spamm_chunk_get_use_linear_tree(chunk_A);
  number_dimensions = *spamm_chunk_get_number_dimensions(chunk_A);

  if(use_linear_tree)
  {
    return spamm_linear_multiply(tolerance, alpha, chunk_A, chunk_B, chunk_C);
  }

  else
  {
    norm_A = spamm_chunk_get_norm(chunk_A);
    norm_B = spamm_chunk_get_norm(chunk_B);

    norm_C = spamm_chunk_get_norm(chunk_C);
    norm2_C = spamm_chunk_get_norm2(chunk_C);

    if(norm_A[0]*norm_B[0] > tolerance)
    {
      matrix_A = spamm_chunk_get_matrix(chunk_A);
      matrix_B = spamm_chunk_get_matrix(chunk_B);
      matrix_C = spamm_chunk_get_matrix(chunk_C);

      N_contiguous = spamm_chunk_get_N_contiguous(chunk_A);

      if(sgemm)
      {
        sgemm("N", "N", &N_contiguous, &N_contiguous, &N_contiguous,
            &alpha_sgemm, matrix_A, &N_contiguous, matrix_B, &N_contiguous,
            &beta, matrix_C, &N_contiguous);
      }

      else
      {
        /* Braindead multiply in nested loops. */
        for(i = 0; i < N_contiguous; i++) {
          for(j = 0; j < N_contiguous; j++) {
            for(k = 0; k < N_contiguous; k++)
            {
              matrix_C[spamm_index_column_major(i, j, N_contiguous, N_contiguous)] += alpha
                *matrix_A[spamm_index_column_major(i, k, N_contiguous, N_contiguous)]
                *matrix_B[spamm_index_column_major(k, j, N_contiguous, N_contiguous)];
            }
          }
        }
      }

      /* Update norm on C. */
      norm2_C[0] = 0;
      for(i = 0; i < ipow(N_contiguous, number_dimensions); i++)
      {
        norm2_C[0] += ipow(matrix_C[i], 2);
      }
      norm_C[0] = sqrt(norm2_C[0]);
    }

    return norm_C[0];
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
  unsigned int *new_N_lower;
  unsigned int *new_N_upper;

  short i, j, k;

  if(node_A == NULL || node_B == NULL) { return; }

  /* We have to allocate a new C block a tier up. */
  if(*node_C == NULL)
  {
    *node_C = spamm_recursive_new_node();
  }

  if(tier == chunk_tier)
  {
    (*node_C)->norm2 = spamm_chunk_multiply(tolerance, alpha,
        node_A->tree.chunk, node_B->tree.chunk, (*node_C)->tree.chunk, sgemm);
    (*node_C)->norm = sqrt((*node_C)->norm2);
  }

  else
  {
    if((*node_C)->tree.child == NULL)
    {
      (*node_C)->tree.child = calloc(ipow(2, number_dimensions_C), sizeof(struct spamm_recursive_node_t*));
    }

    new_N_lower = calloc(number_dimensions_C, sizeof(unsigned int));
    new_N_upper = calloc(number_dimensions_C, sizeof(unsigned int));

    if(number_dimensions_A == 2 &&
        number_dimensions_B == 2 &&
        number_dimensions_C == 2)
    {
      for(i = 0; i < 2; i++) {
        for(j = 0; j < 2; j++) {
          for(k = 0; k < 2; k++)
          {
            new_N_lower[0] = N_lower[0]+(N_upper[0]-N_lower[0])/2*i;
            new_N_upper[0] = N_lower[0]+(N_upper[0]-N_lower[0])/2*(i+1);
            new_N_lower[1] = N_lower[1]+(N_upper[1]-N_lower[1])/2*j;
            new_N_upper[1] = N_lower[1]+(N_upper[1]-N_lower[1])/2*(j+1);

            if(node_A->norm*node_B->norm > tolerance)
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

  if(A->chunk_tier == 0)
  {
    spamm_chunk_multiply_scalar(beta, C->tree.chunk);
    spamm_chunk_multiply(tolerance, alpha, A->tree.chunk, B->tree.chunk,
        C->tree.chunk, sgemm);
  }

  else
  {
    spamm_recursive_multiply_scalar(beta, C->tree.recursive_tree,
        C->number_dimensions, 0, C->chunk_tier, C->use_linear_tree);

    N_lower = calloc(C->number_dimensions, sizeof(unsigned int));
    N_upper = calloc(C->number_dimensions, sizeof(unsigned int));

    for(dim = 0; dim < A->number_dimensions; dim++)
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
