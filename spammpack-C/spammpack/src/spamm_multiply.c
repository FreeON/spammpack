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

/** @private List of index (key) space of kernel tier.
 */
struct spamm_multiply_index_list_t
{
  /** The number of elements in the index array. */
  unsigned int size;

  /** An array that contains a list of linear 2D matrix indices. */
  struct spamm_list_t *index;

  /** An array of pointers to the matrix tree nodes at the kernel tier. */
  struct spamm_hashed_data_t **data;
};

/** @private k lookup table to speed up loop over indices. */
struct spamm_multiply_k_lookup_t
{
  /** The number of elements in the index array. */
  unsigned int size;

  /** An array of index values pointing into the linear index array for
   * different k values. */
  unsigned int *index;
};

/** @private Multiply a matrix node by a scalar.
 *
 * @param index The linear index of that matrix node.
 * @param value A pointer to the spamm_hashed_data_t matrix node.
 * @param user_data The scalar \f$\beta\f$ that multiplies the matrix.
 */
void
spamm_multiply_beta_block (unsigned int index, void *value, void *user_data)
{
  struct spamm_hashed_data_t *data = value;
  float *beta = user_data;
  unsigned int i;

#ifdef HAVE_SSE_DISABLED
  short j;
  __m128 xmm, xmm_beta;

  xmm_beta = _mm_load_ps1(beta);
  for (i = 0; i < SPAMM_N_KERNEL*SPAMM_N_KERNEL; i += 4)
  {
    xmm = _mm_load_ps(&data->block_dense[i]);
    xmm = _mm_mul_ss(xmm_beta, xmm);
    _mm_store_ps(&data->block_dense[i], xmm);

    for (j = 0; j < 4; j++)
    {
      xmm = _mm_load_ps(&data->block_dense_dilated[4*i+4*j]);
      xmm = _mm_mul_ss(xmm_beta, xmm);
      _mm_store_ps(&data->block_dense_dilated[4*i+4*j], xmm);
    }
  }
#else
  for (i = 0; i < SPAMM_N_KERNEL*SPAMM_N_KERNEL; i++)
  {
    /* We only multiply data in block_dense, not the dilated, nor the
     * transpose data (if it exists). This obviously introduces incosistencies
     * into the C matrix, the stream product however, does that right now
     * already anyway. */
    data->block_dense[i] *= (*beta);
  }
#endif
}

void
spamm_hashed_multiply_scalar (const float alpha, struct spamm_hashed_t *A)
{
  struct spamm_hashtable_t *tier_hashtable;

  if (A == NULL) { return; }

  if (alpha != 1.0)
  {
    tier_hashtable = A->tier_hashtable[A->kernel_tier-A->tier];
    spamm_hashtable_foreach(tier_hashtable, spamm_multiply_beta_block, (void*) &alpha);
  }
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

/** @private Swap 2 multiply stream elements.
 *
 * @param a_stream The first stream element.
 * @param b_stream The second stream element.
 * @param a The first linear index of the C matrix.
 * @param b The second linear index of the C matrix.
 */
void
spamm_multiply_sort_stream_swap (struct spamm_multiply_stream_t *a_stream,
    struct spamm_multiply_stream_t *b_stream, unsigned int *a, unsigned int *b)
{
  struct spamm_hashed_data_t *temp_node;
  unsigned int temp;

  temp_node = a_stream->A;
  a_stream->A= b_stream->A;
  b_stream->A= temp_node;

  temp_node = a_stream->B;
  a_stream->B= b_stream->B;
  b_stream->B= temp_node;

  temp_node = a_stream->C;
  a_stream->C= b_stream->C;
  b_stream->C= temp_node;

  temp = *a;
  *a = *b;
  *b = temp;
}

/** @private Sort the multiply stream according to a linear 2D index. This is
 * used to sort the stream according to the linear index of the C blocks to
 * help avoid excessive hash table lookups in associating C blocks to the
 * stream. The sub-list sorted is given by the indices left and right, such
 * that [left, right], i.e. right is inclusive in the array.
 *
 * @param left The left index of the sub-list to be sorted.
 * @param right The right index of the sub-list to be sorted.
 * @param multiply_stream The multiply stream.
 * @param C_block_stream_index The array of linear matrix indices of the C
 * blocks.
 */
void
spamm_multiply_sort_stream (const unsigned int left,
    const unsigned int right,
    struct spamm_multiply_stream_t *multiply_stream,
    unsigned int *C_block_stream_index)
{
  unsigned int i;
  unsigned int pivot;
  unsigned int new_pivot;
  unsigned int pivot_value;

  if (right > left)
  {
    /* Select pivot value. */
    pivot = left+(right-left)/2;

    /* Partition. */
    pivot_value = C_block_stream_index[pivot];

    /* Move pivot to the end. */
    spamm_multiply_sort_stream_swap(&multiply_stream[pivot], &multiply_stream[right],
        &C_block_stream_index[pivot], &C_block_stream_index[right]);

    /* Find new pivot. */
    new_pivot = left;

    for (i = left; i < right; i++)
    {
      if (C_block_stream_index[i] <= pivot_value)
      {
        spamm_multiply_sort_stream_swap(&multiply_stream[i], &multiply_stream[new_pivot], &C_block_stream_index[i], &C_block_stream_index[new_pivot]);
        new_pivot++;
      }
    }
    spamm_multiply_sort_stream_swap(&multiply_stream[new_pivot], &multiply_stream[right], &C_block_stream_index[new_pivot], &C_block_stream_index[right]);

    /* Recurse. */
    spamm_multiply_sort_stream(left, new_pivot-1, multiply_stream, C_block_stream_index);
    spamm_multiply_sort_stream(new_pivot+1, right, multiply_stream, C_block_stream_index);
  }
}

/** Sort the C index array.
 *
 * @param array The array to sort.
 * @param stream The multiply stream.
 * @param length The length of the array and the stream.
 */
void
spamm_multiply_C_index_sort (unsigned int *array,
    struct spamm_multiply_stream_t *stream,
    const unsigned int length)
{
  unsigned int i, j, j_next, i_left, i_right;
  unsigned int sub_current, sub_next;
  unsigned int sub_length;
  unsigned int *sublist;
  unsigned int *scratch_index;
  struct spamm_multiply_stream_t *scratch_stream;

  /* The array is trivially already sorted. */
  if (length <= 1) { return; }

  /* Create index array for sublists. This array is length+1 since we
   * terminate the array by a value of length. */
  sublist = (unsigned int*) malloc(sizeof(unsigned int)*2*(length+1));

  /* Allocate scratch space. */
  scratch_index = (unsigned int*) calloc(length, sizeof(unsigned int));
  scratch_stream = (struct spamm_multiply_stream_t*) calloc(length, sizeof(struct spamm_multiply_stream_t));

  /* Break the original list into at most N pieces, i.e. single element
   * sublists. If adajacent list elements are already in the right order, we
   * put them into the same sublist. */
  sublist[0] = 0;
  for (i = 1, j = 1; i < length; i++)
  {
    if (array[i-1] > array[i])
    {
      /* The 2 elements are in incorrect order. Start a new sublist. */
      sublist[j++] = i;
    }
  }
  sublist[j++] = i;
  sub_length = j;

  /* Loop over the list, merging neighboring sublists until everthying is
   * sorted. */
  sub_current = 0;
  sub_next = length+1;
  while (sublist[sub_current+1] < length)
  {
    for (j = 0, j_next = 0; j < sub_length-2; j += 2)
    {
      /* Merge 2 adjacent sublists. */
      for (i = sublist[sub_current+j], i_left = sublist[sub_current+j], i_right = sublist[sub_current+j+1];
          i < sublist[sub_current+j+2];
          i++)
      {
        if (i_left < sublist[sub_current+j+1] && i_right < sublist[sub_current+j+2])
        {
          if (array[i_left] <= array[i_right])
          {
            scratch_index[i] = array[i_left];
            scratch_stream[i].A = stream[i_left].A;
            scratch_stream[i].B = stream[i_left].B;
            scratch_stream[i].C = stream[i_left].C;
            i_left++;
          }

          else
          {
            scratch_index[i] = array[i_right];
            scratch_stream[i].A = stream[i_right].A;
            scratch_stream[i].B = stream[i_right].B;
            scratch_stream[i].C = stream[i_right].C;
            i_right++;
          }
        }

        else if (i_left < sublist[sub_current+j+1])
        {
          scratch_index[i] = array[i_left];
          scratch_stream[i].A = stream[i_left].A;
          scratch_stream[i].B = stream[i_left].B;
          scratch_stream[i].C = stream[i_left].C;
          i_left++;
        }

        else
        {
          scratch_index[i] = array[i_right];
          scratch_stream[i].A = stream[i_right].A;
          scratch_stream[i].B = stream[i_right].B;
          scratch_stream[i].C = stream[i_right].C;
          i_right++;
        }
      }

      /* Copy the merged list back. */
      for (i = sublist[sub_current+j]; i < sublist[sub_current+j+2]; i++)
      {
        array[i] = scratch_index[i];
        stream[i].A = scratch_stream[i].A;
        stream[i].B = scratch_stream[i].B;
        stream[i].C = scratch_stream[i].C;
      }

      /* Remove division between the sublists just merged. */
      sublist[sub_next+j_next] = sublist[sub_current+j];
      sublist[sub_next+j_next+1] = sublist[sub_current+j+2];
      j_next++;
    }

    /* Add remaining sublist divisions. */
    while (j < sub_length)
    {
      sublist[sub_next+j_next] = sublist[sub_current+j];
      j++;
      j_next++;
    }
    sub_length = j_next;

    /* Switch sublists. */
    if (sub_current == 0) { sub_current = length+1; }
    else                  { sub_current = 0; }
    if (sub_next == 0) { sub_next = length+1; }
    else               { sub_next = 0; }
  }

  /* Free memory. */
  free(sublist);
  free(scratch_index);
  free(scratch_stream);
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

  unsigned int *stream_A;
  unsigned int *stream_B;
  unsigned int *stream_C;

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
  stream_A = malloc(ipow(index_length, 2)*sizeof(unsigned int));
  stream_B = malloc(ipow(index_length, 2)*sizeof(unsigned int));
  stream_C = malloc(ipow(index_length, 2)*sizeof(unsigned int));

  for (i = 0, stream_index = 0; i < index_length; i++) {
    for (j = 0; j < index_length; j++)
    {
      if (norm_A[index_A[i]]*norm_B[index_B[j]] > tolerance)
      {
        stream_A[stream_index] = index_A[i];
        stream_B[stream_index] = index_B[i];
        stream_C[stream_index] = (index_A[i] & MASK_2D_I) | (index_B[j] & MASK_2D_J);
        stream_index++;
      }
    }
  }
  printf("[multiply] Added %u block products to stream\n", stream_index);

  /* Run kernel. */
  norm_A = spamm_chunk_get_tier_norm(*spamm_chunk_get_number_tiers(chunk_A)+SPAMM_KERNEL_DEPTH, chunk_A);
  norm_B = spamm_chunk_get_tier_norm(*spamm_chunk_get_number_tiers(chunk_B)+SPAMM_KERNEL_DEPTH, chunk_B);
  for (i = 0; i < stream_index; i++)
  {
  }

  /* Free memory. */
  free(stream_A);
  free(stream_B);
  free(stream_C);
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
 * @param depth The depth of the matrix tree.
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
    const unsigned int N_block,
    const unsigned int depth,
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
        chunk_tier, depth, N_block, 0, 0, 0, timer, sgemm, kernel);
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
                  number_dimensions_A, number_dimensions_B, number_dimensions_C,
                  N, new_N_lower, new_N_upper, tier+1, chunk_tier, N_block,
                  depth, use_linear_tree, number_products);
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
        C->tree.chunk, 0, C->chunk_tier, C->depth, C->N_block, 0, 0, 0,
        timer, sgemm, kernel);
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
        A->N_block, A->depth, A->use_linear_tree, number_products);

    free(N_lower);
    free(N_upper);
  }
}
