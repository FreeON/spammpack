/** @file */

#include "config.h"
#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef SPAMM_MULTIPLY_DEBUG
#include <string.h>
#endif

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

/** Multiply a chunk with a scalar.
 *
 * @param alpha The factor.
 * @param chunk The chunk.
 * @param flop The flop count.
 *
 * @return The new squared norm.
 */
spamm_norm_t
spamm_chunk_multiply_scalar (const float alpha,
    spamm_chunk_t *chunk,
    double *const flop)
{
  unsigned int i;
  unsigned int number_norms;
  unsigned int number_tiers;
  unsigned int N_contiguous;
  unsigned int number_dimensions;

  float *A;
  float *A_dilated;

  spamm_norm_t *norm;
  spamm_norm_t *norm2;

  if(chunk == NULL) { return 0.0; }

  number_dimensions = *spamm_chunk_get_number_dimensions(chunk);
  number_tiers = *spamm_chunk_get_number_tiers(chunk);
  N_contiguous = spamm_chunk_get_N_contiguous(chunk);
  A = spamm_chunk_get_matrix(chunk);
  A_dilated = spamm_chunk_get_matrix_dilated(chunk);

  for(i = 0; i < ipow(N_contiguous, number_dimensions); i++)
  {
    A[i] *= alpha;

    A_dilated[4*i+0] = A[i];
    A_dilated[4*i+1] = A[i];
    A_dilated[4*i+2] = A[i];
    A_dilated[4*i+3] = A[i];
  }

  /* Update the flop count. */
  *flop += ipow(N_contiguous, number_dimensions);

  norm = spamm_chunk_get_norm(chunk);
  norm2 = spamm_chunk_get_norm2(chunk);
  number_norms = spamm_chunk_get_total_number_norms(number_tiers, number_dimensions);
  for(i = 0; i < number_norms; i++)
  {
    norm[i] *= alpha;
    norm2[i] *= alpha*alpha;
  }

  return norm2[0];
}

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
    const short use_linear_tree,
    double *const flop)
{
  unsigned int i;

  if(A == NULL) { return; }

  if(tier == chunk_tier)
  {
    A->norm2 = spamm_chunk_multiply_scalar(alpha, A->tree.chunk, flop);
    A->norm = sqrt(A->norm2);
  }

  else
  {
    if(A->tree.child != NULL)
    {
      for(i = 0; i < ipow(2, number_dimensions); i++)
      {
        spamm_recursive_multiply_scalar(alpha, A->tree.child[i],
            number_dimensions, tier+1, chunk_tier, use_linear_tree, flop);
      }
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
 * @param flop Number of floating point operations.
 *
 * @return The square of the Frobenius norm of the chunk.
 */
double
spamm_linear_multiply (const spamm_norm_t tolerance,
    const float alpha,
    const spamm_chunk_t *const chunk_A,
    const spamm_chunk_t *const chunk_B,
    spamm_chunk_t *const chunk_C,
    double *const flop)
{
  unsigned int *index_A;
  unsigned int *index_B;

  unsigned int *stream;

  unsigned int N_contiguous;
  unsigned int index_length;

  unsigned int i, j;
  unsigned int i_block, j_block;
  unsigned int j_A, j_B;
  unsigned int stream_index;
  unsigned int i_stream;

  int tier;

  unsigned int norm_offset_C;
  unsigned int offset_C;

  spamm_norm_t *norm_A;
  spamm_norm_t *norm_B;

  spamm_norm_t *norm_C;
  spamm_norm_t *norm2_C;
  spamm_norm_t *norm_C_next;
  spamm_norm_t *norm2_C_next;

  float *matrix_C;
  float *matrix_dilated_C;

#if !defined(RUN_ASSEMBLY_KERNEL)
  float *matrix_A;
  float *matrix_B;

  unsigned int norm_offset_A, norm_offset_B;
  unsigned int offset_A, offset_B;
#endif

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
  SPAMM_WARN("potentially %u product(s)\n", ipow(index_length, 3));
#endif

  /* Convolute. */
  if((stream = malloc(ipow(index_length, 3)*3*sizeof(unsigned int))) == NULL)
  {
    SPAMM_FATAL("can not allocate stream\n");
  }

#ifdef SPAMM_MULTIPLY_DEBUG
  SPAMM_WARN("allocated stream (%p) with %i triple elements\n", stream, ipow(index_length, 3));
#endif

  for(i = 0, stream_index = 0; i < index_length; i++)
  {
    for(j_A = i*index_length; j_A < (i+1)*index_length; j_A++) {
      for(j_B = i*index_length; j_B < (i+1)*index_length; j_B++)
      {
        if(norm_A[index_A[j_A]]*norm_B[index_B[j_B]] > tolerance)
        {
          stream[3*stream_index+0] = index_A[j_A];
          stream[3*stream_index+1] = index_B[j_B];
          stream[3*stream_index+2] = (index_A[j_A] & MASK_2D_I) | (index_B[j_B] & MASK_2D_J);
#ifdef SPAMM_MULTIPLY_DEBUG
          SPAMM_WARN("comparing norms: %e (norm_A[%u]) * %e (norm_B[%u]) = %e -> adding stream[%u] = { %u, %u, %u }\n",
              norm_A[index_A[j_A]],
              index_A[j_A],
              norm_B[index_B[j_B]],
              index_B[j_B],
              norm_A[index_A[j_A]]*norm_B[index_B[j_B]],
              stream_index,
              stream[3*stream_index+0],
              stream[3*stream_index+1],
              stream[3*stream_index+2]);
#endif
          stream_index++;
        }

        else
        {
#ifdef SPAMM_MULTIPLY_DEBUG
          SPAMM_WARN("comparing norms: %e (norm_A[%u]) * %e (norm_B[%u]) = %e -> skipping\n",
              norm_A[index_A[j_A]],
              index_A[j_A],
              norm_B[index_B[j_B]],
              index_B[j_B],
              norm_A[index_A[j_A]]*norm_B[index_B[j_B]],
              stream);
#endif
          break;
        }
      }
    }
  }

  /* Free memory. */
  free(index_A);
  free(index_B);

#ifdef SPAMM_MULTIPLY_DEBUG
  SPAMM_WARN("Added %u (out of %u possible) block products to stream\n", stream_index, ipow(index_length, 3));
  SPAMM_WARN("stream = %p\n", stream);
#endif

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

#ifdef SPAMM_MULTIPLY_DEBUG
                printf("(%u,%u) -> row_maj = %u\n", i_block, j_block, i_block*SPAMM_N_BLOCK+j_block);
                fflush(stdout);
#endif
              }
            }

            /* Update flop count. */
            *flop += 2*SPAMM_N_BLOCK*SPAMM_N_BLOCK*SPAMM_N_BLOCK;
          }
        }
      }
    }
  }
#endif

  /* Free memory. */
#ifdef SPAMM_MULTIPLY_DEBUG
  SPAMM_WARN("freeing stream = %p\n", stream);
#endif
  free(stream);

  /* Update norms. */
  spamm_chunk_fix(chunk_C);

  return norm2_C_next[0];
}

/** Multiply two SpAMM chunks.
 *
 * @param tolerance The SpAMM tolerance.
 * @param alpha The scalar alpha.
 * @param chunk_A Chunk A.
 * @param chunk_B Chunk B.
 * @param chunk_C Chunk C.
 * @param sgemm The sgemm() function to call.
 * @param flop The flop count.
 *
 * @return The square of the norm.
 */
double
spamm_chunk_multiply (const spamm_norm_t tolerance,
    const float alpha,
    const spamm_chunk_t *const chunk_A,
    const spamm_chunk_t *const chunk_B,
    spamm_chunk_t *const chunk_C,
    sgemm_func sgemm,
    double *const flop)
{
  unsigned int i, j, k;

  spamm_norm_t *norm_A;
  spamm_norm_t *norm_B;

  spamm_norm_t *norm_C;
  spamm_norm_t *norm2_C;

  float alpha_sgemm = alpha;
  float beta = 1.0;

  int N_contiguous;

  float *matrix_A;
  float *matrix_B;
  float *matrix_C;

  unsigned int number_dimensions;

  assert(chunk_C != NULL);

  if(chunk_A == NULL || chunk_B == NULL) { return 0.0; }

  number_dimensions = *spamm_chunk_get_number_dimensions(chunk_A);

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

    /* Update flop count. */
    *flop += 2*N_contiguous*N_contiguous*N_contiguous;

    /* Update norm on C. */
    for(i = 0, norm2_C[0] = 0; i < ipow(N_contiguous, number_dimensions); i++)
    {
      norm2_C[0] += matrix_C[i]*matrix_C[i];
    }
    norm_C[0] = sqrt(norm2_C[0]);
  }

  return norm2_C[0];
}

/** Multiply two matrices, i.e. \f$ C = \alpha A \times B + \beta C\f$.
 *
 * @param tolerance The SpAMM tolerance of this product.
 * @param alpha The paramater \f$\alpha\f$.
 * @param A The matrix \f$A\f$.
 * @param B The matrix \f$B\f$.
 * @param beta The paramater \f$\beta\f$.
 * @param C The matrix \f$C\f$.
 * @param sgemm The external sgemm function to use.
 * @param tier The tier.
 * @param chunk_tier The contiguous tier.
 * @param use_linear_tree If set to 1 then we will switch to a linear tree at
 * chunk_tier.
 * @param flop The number of floating point operations.
 */
void
spamm_recursive_multiply (const spamm_norm_t tolerance,
    const float alpha,
    struct spamm_recursive_node_t *node_A,
    struct spamm_recursive_node_t *node_B,
    struct spamm_recursive_node_t *node_C,
    sgemm_func sgemm,
    const unsigned int number_dimensions_A,
    const unsigned int number_dimensions_B,
    const unsigned int number_dimensions_C,
    const unsigned int *const N,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    const unsigned int tier,
    const unsigned int chunk_tier,
    const short use_linear_tree,
    double *const flop)
{
  unsigned int *new_N_lower = NULL;
  unsigned int *new_N_upper = NULL;

  short i, j, k;

#ifdef SPAMM_MULTIPLY_DEBUG
  char temp[100];
  char output_buffer[5000];

  snprintf(output_buffer, 2000-1, "tier = %i, node_A = %p, node_B = %p, node_C = %p, C", tier, node_A, node_B, node_C);
  for(i = 0; i < number_dimensions_C; i++)
  {
    snprintf(temp, 100-1, "[%i->%i]", N_lower[i], N_upper[i]);
    strcat(output_buffer, temp);
  }
  strcat(output_buffer, "\n");
  SPAMM_WARN(output_buffer);

  SPAMM_WARN("&tolerance = %p, node_A = %p, &node_A = %p, node_B = %p, &node_B = %p, node_C = %p, &node_C = %p\n",
      &tolerance, node_A, &node_A, node_B, &node_B, node_C, &node_C);
#endif

  /* node_C has to be set, otherwise we won't have a lock. */
  assert(node_C != NULL);

  if(node_A == NULL || node_B == NULL)
  {
    return;
  }

  if(tier == chunk_tier)
  {
    if(node_A->tree.chunk == NULL || node_B->tree.chunk == NULL)
    {
      return;
    }

#ifdef _OPENMP
    omp_set_lock(&node_C->lock);
#endif

    if(node_C->tree.chunk == NULL)
    {
      node_C->tree.chunk = spamm_new_chunk(number_dimensions_C, use_linear_tree, N, N_lower, N_upper);
    }

    if(use_linear_tree)
    {
      node_C->norm2 = spamm_linear_multiply(tolerance, alpha,
          node_A->tree.chunk, node_B->tree.chunk, node_C->tree.chunk, flop);
    }

    else
    {
      node_C->norm2 = spamm_chunk_multiply(tolerance, alpha,
          node_A->tree.chunk, node_B->tree.chunk, node_C->tree.chunk, sgemm, flop);
    }

    node_C->norm = sqrt(node_C->norm2);

#ifdef _OPENMP
    omp_unset_lock(&node_C->lock);
#endif
  }

  else
  {
    if(node_A->tree.child == NULL || node_B->tree.child == NULL)
    {
      return;
    }

    if(number_dimensions_A == 2 &&
       number_dimensions_B == 2 &&
       number_dimensions_C == 2)
    {
      for(k = 0; k < 2; k++) {
        for(i = 0; i < 2; i++) {
          for(j = 0; j < 2; j++)
          {
            if(node_A->tree.child[i+2*k] != NULL && node_B->tree.child[k+2*j] != NULL)
            {
              if(node_A->tree.child[i+2*k]->norm*node_B->tree.child[k+2*j]->norm > tolerance)
              {
#ifdef _OPENMP
                omp_set_lock(&node_C->lock);
#endif
                if(node_C->tree.child == NULL)
                {
                  node_C->tree.child = calloc(ipow(2, number_dimensions_C), sizeof(struct spamm_recursive_node_t*));
                }

                if(node_C->tree.child[i+2*j] == NULL)
                {
                  node_C->tree.child[i+2*j] = spamm_recursive_new_node();
                }
#ifdef _OPENMP
                omp_unset_lock(&node_C->lock);
#endif

#ifdef SPAMM_MULTIPLY_DEBUG
                SPAMM_WARN("descending: C[%i][%i] (%p) <- A[%i][%i]*B[%i][%i] (%e)\n", i, j, node_C->tree.child[i+2*j],
                    i, k, k, j, node_A->tree.child[i+2*k]->norm*node_B->tree.child[k+2*j]->norm);
#endif

#pragma omp task untied
                {
                  new_N_lower = calloc(number_dimensions_C, sizeof(unsigned int));
                  new_N_upper = calloc(number_dimensions_C, sizeof(unsigned int));

                  new_N_lower[0] = N_lower[0]+(N_upper[0]-N_lower[0])/2*i;
                  new_N_upper[0] = N_lower[0]+(N_upper[0]-N_lower[0])/2*(i+1);
                  new_N_lower[1] = N_lower[1]+(N_upper[1]-N_lower[1])/2*j;
                  new_N_upper[1] = N_lower[1]+(N_upper[1]-N_lower[1])/2*(j+1);

#ifdef SPAMM_MULTIPLY_DEBUG
                  SPAMM_WARN("created thread: node_C = %p, node_C->tree.child[%i] C[%i->%i][%i->%i] = %p\n",
                      node_C, i+2*j,
                      new_N_lower[0], new_N_upper[0], new_N_lower[1], new_N_upper[1],
                      node_C->tree.child[i+2*j]);
#endif

                  spamm_recursive_multiply(tolerance, alpha,
                      node_A->tree.child[i+2*k], node_B->tree.child[k+2*j],
                      node_C->tree.child[i+2*j], sgemm, number_dimensions_A,
                      number_dimensions_B, number_dimensions_C, N, new_N_lower,
                      new_N_upper, tier+1, chunk_tier, use_linear_tree,
                      flop);

                  free(new_N_lower);
                  free(new_N_upper);
                }
              }

#ifdef SPAMM_MULTIPLY_DEBUG
              else
              {
                SPAMM_WARN("skipping product below tolerance: C[%i][%i] <- A[%i][%i]*B[%i][%i] (%e)\n", i, j, i, k, k, j,
                    node_A->tree.child[i+2*k]->norm*node_B->tree.child[k+2*j]->norm);
              }
#endif
            }
          }
        }
      }
#pragma omp taskwait

      /* Fix up norms. */
      if(node_C->tree.child != NULL)
      {
        for(i = 0, node_C->norm2 = 0; i < ipow(2, number_dimensions_C); i++)
        {
          if(node_C->tree.child[i] != NULL)
          {
            node_C->norm2 += node_C->tree.child[i]->norm2;
          }
        }
        node_C->norm = sqrt(node_C->norm2);
      }
    }

    else
    {
      SPAMM_FATAL("not implemented\n");
    }
  }
}

/** Multiply a matrix by a scalar, @f$ A \leftarrow \alpha A @f$.
 *
 * @param alpha The scalar.
 * @param A The matrix.
 * @param flop The flop count.
 */
void
spamm_multiply_scalar (const float alpha,
    struct spamm_matrix_t *const A,
    double *const flop)
{
  spamm_recursive_multiply_scalar(alpha, A->recursive_tree,
      A->number_dimensions, 0, A->chunk_tier, A->use_linear_tree, flop);
}

/** Multiply two matrices, i.e. \f$ C = \alpha A \times B + \beta C\f$.
 *
 * @param tolerance The SpAMM tolerance of this product.
 * @param alpha The paramater \f$\alpha\f$.
 * @param A The matrix \f$A\f$.
 * @param B The matrix \f$B\f$.
 * @param beta The paramater \f$\beta\f$.
 * @param C The matrix \f$C\f$.
 * @param sgemm The external sgemm function to use.
 * @param flop The number of floating point operations.
 */
void
spamm_multiply (const spamm_norm_t tolerance,
    const float alpha,
    const struct spamm_matrix_t *const A,
    const struct spamm_matrix_t *const B,
    const float beta,
    struct spamm_matrix_t *const C,
    sgemm_func sgemm,
    double *const flop)
{
  int dim;
  unsigned int *N_lower;
  unsigned int *N_upper;

  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);
  assert(flop != NULL);

  spamm_recursive_multiply_scalar(beta, C->recursive_tree,
      C->number_dimensions, 0, C->chunk_tier, C->use_linear_tree, flop);

  if(alpha != 0.0)
  {
    /* Allocate a new tree node before we enter the recursive portion of the
     * multiply. Note that we need to do this in the serial section of the
     * code to prevent a race condition.
     */
    if(C->recursive_tree == NULL)
    {
      C->recursive_tree = spamm_recursive_new_node();
    }

#pragma omp parallel
    {
#pragma omp single
      {
#pragma omp task untied
        {
          /* Allocate and set within task region. Otherwise we will end up
           * with wrong N_{upper,lower} values in the
           * spamm_recursive_multiply() call. */
          N_lower = calloc(C->number_dimensions, sizeof(unsigned int));
          N_upper = calloc(C->number_dimensions, sizeof(unsigned int));

          for(dim = 0; dim < A->number_dimensions; dim++)
          {
            N_upper[dim] = A->N_padded;
          }

          spamm_recursive_multiply(tolerance, alpha, A->recursive_tree,
              B->recursive_tree, C->recursive_tree, sgemm, A->number_dimensions,
              B->number_dimensions, C->number_dimensions, A->N, N_lower, N_upper, 0,
              A->chunk_tier, A->use_linear_tree, flop);

          free(N_lower);
          free(N_upper);
        }
      }
    }

    /* Prune tree. */
    spamm_prune(C);
  }
}
