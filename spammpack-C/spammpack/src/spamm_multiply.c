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

  if (number_dimensions == 2 && use_linear_tree && tier == chunk_tier)
  {
    spamm_hashed_multiply_scalar(alpha, A->tree.hashed_tree);
  }

  else if (tier == chunk_tier)
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
 */
void
spamm_linear_multiply (const float tolerance,
    const float alpha,
    struct spamm_hashed_t *A,
    struct spamm_hashed_t *B,
    spamm_chunk_t *chunk_A,
    spamm_chunk_t *chunk_B,
    const float beta,
    struct spamm_hashed_t *C,
    spamm_chunk_t *chunk_C,
    struct spamm_timer_t *timer,
    const enum spamm_kernel_t kernel)
{
  struct spamm_hashtable_t *A_tier_hashtable;
  struct spamm_hashtable_t *B_tier_hashtable;
  struct spamm_hashtable_t *C_tier_hashtable;

  struct spamm_multiply_index_list_t A_index;
  struct spamm_multiply_index_list_t B_index;

  unsigned int i, j, k, k_check;
  unsigned int index;

  unsigned int A_k_lookup_index;
  unsigned int B_k_lookup_index;
  unsigned int A_k, B_k;

  unsigned int N_contiguous_A;
  unsigned int N_contiguous_B;
  unsigned int N_contiguous_C;

  float *norm_A;
  float *norm_B;
  float *norm_C;

  float *matrix_A;
  float *matrix_B;
  float *matrix_C;

  struct spamm_multiply_k_lookup_t A_k_lookup;
  struct spamm_multiply_k_lookup_t B_k_lookup;

  struct spamm_multiply_stream_t *multiply_stream;
  unsigned int stream_index;
  unsigned int number_dropped_blocks;

  unsigned int *C_block_stream_index;

#ifdef SPAMM_MULTIPLY_PRODUCT_COUNT
  unsigned int number_products = 0;
#endif

  char timer_info_string[2000];

  char *timer_string;

  unsigned int N_contiguous;

  /* For convenience. */
  N_contiguous = spamm_chunk_get_N_contiguous(chunk_A);

  norm_A = spamm_chunk_get_tier_norm(*spamm_chunk_get_number_tiers(chunk_A), chunk_A);
  norm_B = spamm_chunk_get_tier_norm(*spamm_chunk_get_number_tiers(chunk_B), chunk_B);

  matrix_A = spamm_chunk_get_matrix(chunk_A);
  matrix_B = spamm_chunk_get_matrix(chunk_B);

  /* Print out some information. */
#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("[multiply] alpha = %e, beta = %e, tolerance = %e\n", alpha, beta, tolerance);
#endif

  /* Print out some timer information. */
  spamm_timer_info(timer, timer_info_string, 2000);
#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("[multiply] timer: %s\n", timer_info_string);
#endif

#ifdef SPAMM_MULTIPLY_SORT_INDEX
  /* Sort 2D indices on k, i.e. either on row or column index. */
#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("[multiply] sorting A and B on index... ");
#endif
  spamm_timer_start(timer);

  A_index.index = spamm_list_new(N_contiguous_A/SPAMM_N_KERNEL);
  B_index.index = spamm_list_new(N_contiguous_B/SPAMM_N_KERNEL);

  for (i = 0; i < N_contiguous_A/SPAMM_N_KERNEL; i++)
  {
    spamm_list_set(A_index.index, i, i, norm_A[i]);
    spamm_list_set(B_index.index, i, i, norm_B[i]);
  }

  A_index.size = spamm_list_length(A_index.index);
  B_index.size = spamm_list_length(B_index.index);

  spamm_list_sort_index(A_index.index, spamm_list_compare_index_column);
  spamm_list_sort_index(B_index.index, spamm_list_compare_index_row);

#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("len(A) = %u, len(B) = %u, ", spamm_list_length(A_index.index), spamm_list_length(B_index.index));
#endif

  spamm_timer_stop(timer);
  timer_string = spamm_timer_get_string(timer);
#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("%s timer units\n", timer_string);
#endif
  free(timer_string);
#endif

#ifdef SPAMM_MULTIPLY_K_LOOKUP
  /* Create a lookup table for the start of a particular k index in the sorted
   * arrays. */
#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("[multiply] creating k lookup tables... ");
#endif
  spamm_timer_start(timer);

  A_k_lookup.index = spamm_allocate(sizeof(unsigned int)*(N_contiguous/SPAMM_N_KERNEL+1), 1);
  B_k_lookup.index = spamm_allocate(sizeof(unsigned int)*(N_contiguous/SPAMM_N_KERNEL+1), 1);

  /* The index in A_k_lookup. */
  A_k_lookup.size = 0;

  /* The index in A_index.index. */
  i = 0;

  /* The last k value. Initially place it behind the largest expected k value. */
  k = N_contiguous+1;

  for (i = 0; i < spamm_list_length(A_index.index); i++)
  {
    /* Extract k index from linear index. */
    index = spamm_list_get_index(A_index.index, i);
    spamm_index_2D_to_ij(index, NULL, &k_check);

    if (k != k_check)
    {
      A_k_lookup.index[A_k_lookup.size++] = i;
      k = k_check;
    }
  }

  /* Add terminating entry to lookup list. */
  A_k_lookup.index[A_k_lookup.size++] = spamm_list_length(A_index.index);

  /* Check. */
  if (A_k_lookup.size > N_contiguous/SPAMM_N_KERNEL+1)
  {
    SPAMM_FATAL("k lookup table too long for A, estimated %u elements, but found %u\n", N_contiguous/SPAMM_N_KERNEL+1, A_k_lookup.size);
  }

  /* The index in B_k_lookup. */
  B_k_lookup.size = 0;

  /* The index in B_index_sorted. */
  i = 0;

  /* The last k value. Initially place it behind the largest expected k value. */
  k = N_contiguous+1;

  for (i = 0; i < spamm_list_length(B_index.index); i++)
  {
    /* Extract k index from linear index. */
    index = spamm_list_get_index(B_index.index, i);
    spamm_index_2D_to_ij(index, &k_check, NULL);

    if (k != k_check)
    {
      B_k_lookup.index[B_k_lookup.size++] = i;
      k = k_check;
    }
  }

  /* Add terminating entry to lookup list. */
  B_k_lookup.index[B_k_lookup.size++] = spamm_list_length(B_index.index);

  /* Check. */
  if (B_k_lookup.size > N_contiguous/SPAMM_N_KERNEL+1)
  {
    SPAMM_FATAL("k lookup table too long for B, estimated %u elements, but found %u\n", N_contiguous/SPAMM_N_KERNEL+1, B_k_lookup.size);
  }

#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("len(A_k) = %u, len(B_k) = %u, ", A_k_lookup.size, B_k_lookup.size);
#endif

  spamm_timer_stop(timer);
  timer_string = spamm_timer_get_string(timer);
#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("%s timer units\n", timer_string);
#endif
  free(timer_string);
#endif

#ifdef SPAMM_MULTIPLY_SORT_NORM
  /* Sort 2D indices on norms. */
#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("[multiply] sorting A and B on norms... ");
#endif
  spamm_timer_start(timer);

  for (i = 0; i < A_k_lookup.size-1; i++)
  {
    spamm_list_sort_norm(A_index.index, A_k_lookup.index[i], A_k_lookup.index[i+1]);
  }

  for (i = 0; i < B_k_lookup.size-1; i++)
  {
    spamm_list_sort_norm(B_index.index, B_k_lookup.index[i], B_k_lookup.index[i+1]);
  }

  spamm_timer_stop(timer);
  timer_string = spamm_timer_get_string(timer);
#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("%s timer units\n", timer_string);
#endif
  free(timer_string);
#endif

#ifdef SPAMM_MULTIPLY_CONVOLUTE
  /* Convolute by constructing product 3D index. */
#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("[multiply] convolute... ");
#endif
  spamm_timer_start(timer);

  multiply_stream = spamm_allocate(sizeof(struct spamm_multiply_stream_t)
      *(N_contiguous/SPAMM_N_KERNEL)*(N_contiguous/SPAMM_N_KERNEL)*(N_contiguous/SPAMM_N_KERNEL), 0);
  C_block_stream_index = spamm_allocate(sizeof(unsigned int)
      *(N_contiguous/SPAMM_N_KERNEL)*(N_contiguous/SPAMM_N_KERNEL)*(N_contiguous/SPAMM_N_KERNEL), 0);

  unsigned int *A_index_array = spamm_list_get_index_address(A_index.index);
  unsigned int *B_index_array_original = spamm_list_get_index_address(B_index.index);

  unsigned int *multiply_stream_C_index = spamm_allocate(sizeof(unsigned int)
      *(N_contiguous/SPAMM_N_KERNEL)*(N_contiguous/SPAMM_N_KERNEL)*(N_contiguous/SPAMM_N_KERNEL), 0);

  __m128i A_index_xmm;
  __m128i B_index_xmm;
  __m128i C_index_xmm;

  unsigned int C_index_current[4];

  __m128 A_norm_xmm;
  __m128 B_norm_xmm;

  unsigned int *B_k_lookup_array;
  unsigned int *B_index_array;
  unsigned int *B_index_translation;

  float *A_norm_array;
  float *B_norm_array;

  int norm_product[4];
  __m128 norm_product_xmm;

  short int early_termination = 0;

  unsigned int last_C_index;
  struct spamm_hashed_data_t *last_C;

  __m128 tolerance_xmm;

  __m128i mask_2d_i_xmm;
  __m128i mask_2d_j_xmm;

  stream_index = 0;
  number_dropped_blocks = 0;

  tolerance_xmm = _mm_load_ss(&tolerance);
  tolerance_xmm = _mm_shuffle_ps(tolerance_xmm, tolerance_xmm, 0x00);

  mask_2d_i_xmm = _mm_cvtsi32_si128(MASK_2D_I);
  mask_2d_j_xmm = _mm_cvtsi32_si128(MASK_2D_J);

  mask_2d_i_xmm = (__m128i) _mm_shuffle_ps((__m128) mask_2d_i_xmm, (__m128) mask_2d_i_xmm, 0x00);
  mask_2d_j_xmm = (__m128i) _mm_shuffle_ps((__m128) mask_2d_j_xmm, (__m128) mask_2d_j_xmm, 0x00);

  A_norm_array = spamm_allocate(A_index.size*sizeof(float), 0);
  for (i = 0; i < A_index.size; i++)
  {
    A_norm_array[i] = A_index.data[i]->node_norm;
  }

  /* Copy index and norm arrays. This needs to be done to create properly
   * aligned arrays. */
  B_k_lookup_array = spamm_allocate(B_k_lookup.size*sizeof(unsigned int), 1);
  B_norm_array = spamm_allocate(B_index.size*sizeof(float)+B_k_lookup.size*16, 1);
  B_index_array = spamm_allocate(B_index.size*sizeof(unsigned int)+B_k_lookup.size*16, 1);
  B_index_translation = spamm_allocate(B_index.size*sizeof(unsigned int)+B_k_lookup.size*16, 1);

  for (B_k_lookup_index = 0, j = 0; B_k_lookup_index+1 < B_k_lookup.size; B_k_lookup_index++)
  {
    B_k_lookup_array[B_k_lookup_index] = j;

    /* Copy k-block data. */
    for (i = B_k_lookup.index[B_k_lookup_index]; i < B_k_lookup.index[B_k_lookup_index+1]; i++)
    {
      if (i < B_index.size)
      {
        B_norm_array[j] = B_index.data[i]->node_norm;
        B_index_array[j] = B_index_array_original[i];
        B_index_translation[j] = i;
      }

      j++;
    }

    /* Align up to next 128-bit boundary. */
    j = (j+16/sizeof(unsigned int)) & (-16/sizeof(unsigned int));
  }
  B_k_lookup_array[B_k_lookup_index] = j;

  /* Start convolution. */
  for (A_k_lookup_index = 0, B_k_lookup_index = 0; A_k_lookup_index+1 < A_k_lookup.size && B_k_lookup_index+1 < B_k_lookup.size; )
  {
    /* Match k-indices. */
    A_k = A_index_array[A_k_lookup.index[A_k_lookup_index]] & MASK_2D_J;
    B_k = (B_index_array[B_k_lookup_array[B_k_lookup_index]] & MASK_2D_I) >> 1;

    if (A_k > B_k)
    {
      B_k_lookup_index++;
      continue;
    }

    else if (A_k < B_k)
    {
      A_k_lookup_index++;
      continue;
    }

    /* Loop over k-block in A. */
    for (i = A_k_lookup.index[A_k_lookup_index]; i < A_k_lookup.index[A_k_lookup_index+1]; i++)
    {
      /* Get index of A. */
      A_index_xmm = _mm_cvtsi32_si128(A_index_array[i]); /* movd A_index_array_original[i], xmm */
      A_index_xmm = (__m128i) _mm_shuffle_ps((__m128) A_index_xmm, (__m128) A_index_xmm, 0x00); /* pshufd xmm, $0x0 */

      /* Load norm of A. */
      A_norm_xmm = _mm_load_ss(&A_norm_array[i]); /* movss A_norm_current, xmm */
      A_norm_xmm = _mm_shuffle_ps(A_norm_xmm, A_norm_xmm, 0x00); /* shufps xmm, xmm, $0x00 */

      for (j = B_k_lookup_array[B_k_lookup_index], early_termination = 1; j < B_k_lookup_array[B_k_lookup_index+1]; j += 4)
      {
        /* Get index of B. */
        B_index_xmm = _mm_load_si128((__m128i*) &B_index_array[j]); /* movdqa B_index_array[j], xmm */

        /* Load norm of B. */
        B_norm_xmm = _mm_load_ps(&B_norm_array[j]); /* movaps B_norm_array[j], xmm */

        /* Calculate norm products. */
        norm_product_xmm = _mm_mul_ps(A_norm_xmm, B_norm_xmm);
        norm_product_xmm = _mm_cmpgt_ps(norm_product_xmm, tolerance_xmm);
        _mm_store_si128((__m128i*) &norm_product, (__m128i) norm_product_xmm);

        /* Calculate indices of C. */
        C_index_xmm = (__m128i) (_mm_or_ps(
              _mm_and_ps((__m128) A_index_xmm, (__m128) mask_2d_i_xmm),
              _mm_and_ps((__m128) B_index_xmm, (__m128) mask_2d_j_xmm)));
        _mm_store_si128((__m128i*) &C_index_current, C_index_xmm);

        /* Append to stream or done. */
        if (norm_product[0] != 0)
        {
          multiply_stream[stream_index].A = A_index.data[i];
          multiply_stream[stream_index].B = B_index.data[B_index_translation[j]];
          multiply_stream_C_index[stream_index] = C_index_current[0];
          stream_index++;

          /* Since we have been able to find at least one product that passed
           * the SpAMM condition, we have to consider the next block in A. */
          early_termination = 0;
        }

        else
        {
          number_dropped_blocks += B_k_lookup.index[B_k_lookup_index+1]-(j);
          break;
        }

        if (norm_product[1] != 0)
        {
          multiply_stream[stream_index].A = A_index.data[i];
          multiply_stream[stream_index].B = B_index.data[B_index_translation[j+1]];
          multiply_stream_C_index[stream_index] = C_index_current[1];
          stream_index++;
        }

        else
        {
          number_dropped_blocks += B_k_lookup.index[B_k_lookup_index+1]-(j+1);
          break;
        }

        if (norm_product[2] != 0)
        {
          multiply_stream[stream_index].A = A_index.data[i];
          multiply_stream[stream_index].B = B_index.data[B_index_translation[j+2]];
          multiply_stream_C_index[stream_index] = C_index_current[2];
          stream_index++;
        }

        else
        {
          number_dropped_blocks += B_k_lookup.index[B_k_lookup_index+1]-(j+2);
          break;
        }

        if (norm_product[3] != 0)
        {
          multiply_stream[stream_index].A = A_index.data[i];
          multiply_stream[stream_index].B = B_index.data[B_index_translation[j+3]];
          multiply_stream_C_index[stream_index] = C_index_current[3];
          stream_index++;
        }

        else
        {
          number_dropped_blocks += B_k_lookup.index[B_k_lookup_index+1]-(j+3);
          break;
        }
      }

      /* Check whether we had at least one block in B that made it past the
       * SpAMM condition. */
      if (early_termination == 1)
      {
        number_dropped_blocks += (A_k_lookup.index[A_k_lookup_index+1]-(i+1))*(B_k_lookup.index[B_k_lookup_index+1]-B_k_lookup.index[B_k_lookup_index]);
        break;
      }
    }

    A_k_lookup_index++;
  }

  /* Sort multiply_stream and lookup C blocks. */
  //spamm_multiply_C_index_sort(multiply_stream_C_index, multiply_stream, stream_index);

  last_C_index = multiply_stream_C_index[0];
  last_C = spamm_hashtable_lookup(C_tier_hashtable, last_C_index);
  for (i = 0; i < stream_index; i++)
  {
    if (multiply_stream_C_index[i] != last_C_index)
    {
      last_C_index = multiply_stream_C_index[i];
      last_C = spamm_hashtable_lookup(C_tier_hashtable, last_C_index);
    }

    if (last_C == NULL)
    {
      //printf("[%s:%i] creating new node of C\n", __FILE__, __LINE__);
      last_C = spamm_hashed_new_data(C->kernel_tier, last_C_index, C->layout);
      spamm_hashtable_insert(C_tier_hashtable, last_C_index, last_C);
    }

    /* Assign C block. */
    multiply_stream[i].C = last_C;
  }

  /* Free some memory. */
  free(A_norm_array);
  free(B_k_lookup_array);
  free(B_norm_array);
  free(B_index_array);
  free(B_index_translation);
  free(multiply_stream_C_index);

  /* Check. */
  if (stream_index > (N_contiguous/SPAMM_N_KERNEL)*(N_contiguous/SPAMM_N_KERNEL)*(N_contiguous/SPAMM_N_KERNEL))
  {
    SPAMM_FATAL("multiply stream has too many elements, has %u but is only dimensioned for %u\n", stream_index,
        (N_contiguous/SPAMM_N_KERNEL)*(N_contiguous/SPAMM_N_KERNEL)*(N_contiguous/SPAMM_N_KERNEL));
  }

  spamm_timer_stop(timer);
  timer_string = spamm_timer_get_string(timer);
#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("%s timer units\n", timer_string);
#endif
  free(timer_string);
#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("[multiply] dropped %u blocks, placed %u blocks into stream, total of %u blocks\n", number_dropped_blocks, stream_index, number_dropped_blocks+stream_index);
#endif
#endif

#ifdef SPAMM_MULTIPLY_FREE
  /* Free memory. */
#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("[multiply] free memory... ");
#endif
  spamm_timer_start(timer);

  free(A_k_lookup.index);
  free(B_k_lookup.index);

  free(C_block_stream_index);
  spamm_list_delete(&A_index.index);
  free(A_index.data);
  spamm_list_delete(&B_index.index);
  free(B_index.data);

  spamm_timer_stop(timer);
  timer_string = spamm_timer_get_string(timer);
#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("%s timer units\n", timer_string);
#endif
  free(timer_string);
#endif

#ifdef SPAMM_MULTIPLY_STREAM
  /* Call stream product. */
#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("[multiply] stream multiply ");
#endif
  spamm_timer_start(timer);

to be continued...

  spamm_timer_stop(timer);
  timer_string = spamm_timer_get_string(timer);
#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("%s timer units\n", timer_string);
#endif
  free(timer_string);
#endif

#ifdef SPAMM_MULTIPLY_UPDATE_NORM
#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("[multiply] updating C tree norms... ");
#endif
  spamm_timer_start(timer);

  spamm_hashed_norm_update(C);
  spamm_construct_tree(C);

  spamm_timer_stop(timer);
  timer_string = spamm_timer_get_string(timer);
#ifdef SPAMM_MULTIPLY_PRINT_ALOT
  printf("%s timer units\n", timer_string);
#endif
  free(timer_string);
#endif

#ifdef SPAMM_MULTIPLY_PRODUCT_COUNT
  /* Loop through the stream and determine the number of products. */
  printf("[multiply] counting number of block products... ");
  fflush(stdout);
  for (index = 0; index < stream_index; index++)
  {
    for (i = 0; i < SPAMM_N_KERNEL_BLOCKED; i++) {
      for (j = 0; j < SPAMM_N_KERNEL_BLOCKED; j++) {
        for (k = 0; k < SPAMM_N_KERNEL_BLOCKED; k++)
        {
          if (multiply_stream[index].A->norm[spamm_index_norm(i, k)]*multiply_stream[index].B->norm[spamm_index_norm(k, j)] > tolerance)
          {
            number_products++;
          }
        }
      }
    }
  }
  printf("%u products out of %u possible (%1.2f%%)\n", number_products,
      (int) round(ceil(A->N/(double) SPAMM_N_BLOCK)*ceil(A->N/(double) SPAMM_N_BLOCK)*ceil(A->N/(double) SPAMM_N_BLOCK)),
      (double) number_products/round(ceil(A->N/(double) SPAMM_N_BLOCK)*ceil(A->N/(double) SPAMM_N_BLOCK)*ceil(A->N/(double) SPAMM_N_BLOCK))*100);
#endif

#ifdef SPAMM_MULTIPLY_FINAL_FREE
  /* Free memory. */
  free(multiply_stream);
#endif
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
