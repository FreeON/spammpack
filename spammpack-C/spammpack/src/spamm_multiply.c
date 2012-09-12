/** @file */

#include "config.h"
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

/** @private Multiply a matrix by a scalar.
 *
 * @param beta The scalar \f$\beta\f$ that multiplies the matrix.
 * @param A The matrix.
 */
void
spamm_multiply_beta (const float beta, struct spamm_hashed_t *A)
{
  struct spamm_hashtable_t *tier_hashtable;

  if (beta != 1.0)
  {
    tier_hashtable = A->tier_hashtable[A->kernel_tier];
    spamm_hashtable_foreach(tier_hashtable, spamm_multiply_beta_block, (void*) &beta);
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
spamm_hashed_multiply (const float tolerance,
    const float alpha, struct spamm_hashed_t *A, struct spamm_hashed_t *B,
    const float beta, struct spamm_hashed_t *C,
    struct spamm_timer_t *timer,
    const enum spamm_kernel_t kernel)
{
  struct spamm_hashtable_t *A_tier_hashtable;
  struct spamm_hashtable_t *B_tier_hashtable;
  struct spamm_hashtable_t *C_tier_hashtable;

  struct spamm_multiply_index_list_t A_index;
  struct spamm_multiply_index_list_t B_index;

#if SPAMM_MULTIPLY_CONVOLUTE_IMPLEMENTATION == 1
  struct spamm_hashed_data_t *A_data;
  struct spamm_hashed_data_t *B_data;
  struct spamm_hashed_data_t *C_data;
#endif

  unsigned int i, j, k, k_check;
  unsigned int index;
#if SPAMM_MULTIPLY_CONVOLUTE_IMPLEMENTATION == 1
  unsigned int convolution_index_2D;
#endif
  unsigned int A_k_lookup_index;
  unsigned int B_k_lookup_index;
  unsigned int A_k, B_k;
  unsigned int first_A_k, first_B_k;
  unsigned int next_A_k, next_B_k;

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

  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);

  /* Check some more things. */
  if (A->layout != B->layout || A->layout != C->layout)
  {
    printf("[multiply] inconsistent layout in matrices\n");
    exit(1);
  }

  if (A->layout != spamm_kernel_suggest_layout(kernel))
  {
    printf("[multiply] wrong layout for chosen kernel\n");
    exit(1);
  }

  /* Print out some information. */
  printf("[multiply] alpha = %e, beta = %e, tolerance = %e\n", alpha, beta, tolerance);

  /* Print out some timer information. */
  spamm_timer_info(timer, timer_info_string, 2000);
  printf("[multiply] timer: %s\n", timer_info_string);

#ifdef SPAMM_MULTIPLY_BETA
  /* Multiply C with beta. */
  printf("[multiply] multiplying C with beta... ");
  spamm_timer_start(timer);

  spamm_multiply_beta(beta, C);

  spamm_timer_stop(timer);
  timer_string = spamm_timer_get_string(timer);
  printf("%s timer units\n", timer_string);
  free(timer_string);
#endif

#ifdef SPAMM_MULTIPLY_SORT_INDEX
  /* Sort 2D indices on k, i.e. either on row or column index. */
  printf("[multiply] sorting A and B on index... ");
  spamm_timer_start(timer);

  A_tier_hashtable = A->tier_hashtable[A->kernel_tier];
  B_tier_hashtable = B->tier_hashtable[B->kernel_tier];
  C_tier_hashtable = C->tier_hashtable[C->kernel_tier];

  A_index.index = spamm_hashtable_keys(A_tier_hashtable);
  B_index.index = spamm_hashtable_keys(B_tier_hashtable);

  A_index.size = spamm_list_length(A_index.index);
  B_index.size = spamm_list_length(B_index.index);

  spamm_list_sort_index(A_index.index, spamm_list_compare_index_column);
  spamm_list_sort_index(B_index.index, spamm_list_compare_index_row);

  printf("len(A) = %u, len(B) = %u, ", spamm_list_length(A_index.index), spamm_list_length(B_index.index));

  spamm_timer_stop(timer);
  timer_string = spamm_timer_get_string(timer);
  printf("%s timer units\n", timer_string);
  free(timer_string);
#endif

#ifdef SPAMM_MULTIPLY_K_LOOKUP
  /* Create a lookup table for the start of a particular k index in the sorted
   * arrays. */
  printf("[multiply] creating k lookup tables... ");
  spamm_timer_start(timer);

  A_k_lookup.index = spamm_allocate(sizeof(unsigned int)*(A->N_padded/SPAMM_N_KERNEL+1));
  B_k_lookup.index = spamm_allocate(sizeof(unsigned int)*(B->N_padded/SPAMM_N_KERNEL+1));

  /* The index in A_k_lookup. */
  A_k_lookup.size = 0;

  /* The index in A_index.index. */
  i = 0;

  /* The last k value. Initially place it behind the largest expected k value. */
  k = A->N+1;

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
  if (A_k_lookup.size > A->N_padded/SPAMM_N_KERNEL+1)
  {
    printf("k lookup table too long for A, estimated %u elemens, but found %u\n", A->N_padded/SPAMM_N_KERNEL+1, A_k_lookup.size);
    exit(1);
  }

  /* The index in B_k_lookup. */
  B_k_lookup.size = 0;

  /* The index in B_index_sorted. */
  i = 0;

  /* The last k value. Initially place it behind the largest expected k value. */
  k = B->M+1;

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
  if (B_k_lookup.size > B->N_padded/SPAMM_N_KERNEL+1)
  {
    printf("k lookup table too long for B, estimated %u elemens, but found %u\n", B->N_padded/SPAMM_N_KERNEL+1, B_k_lookup.size);
    exit(1);
  }

  printf("len(A_k) = %u, len(B_k) = %u, ", A_k_lookup.size, B_k_lookup.size);

  spamm_timer_stop(timer);
  timer_string = spamm_timer_get_string(timer);
  printf("%s timer units\n", timer_string);
  free(timer_string);
#endif

#ifdef SPAMM_MULTIPLY_SORT_NORM
  /* Sort 2D indices on norms. */
  printf("[multiply] sorting A and B on norms... ");
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
  printf("%s timer units\n", timer_string);
  free(timer_string);
#endif

#ifdef SPAMM_MULTIPLY_COPY_INDICES
  /* Copy sorted indices to array for quick access. */
  printf("[multiply] copying indices to array and referencing dense blocks... ");
  spamm_timer_start(timer);

  A_index.data = spamm_allocate(sizeof(struct spamm_hashed_data_t*)*spamm_list_length(A_index.index));
  B_index.data = spamm_allocate(sizeof(struct spamm_hashed_data_t*)*spamm_list_length(B_index.index));

  printf("len(A_index) = %u, len(B_index) = %u, ", A_index.size, B_index.size);

  for (i = 0; i < A_index.size; i++)
  {
    A_index.data[i] = spamm_hashtable_lookup(A_tier_hashtable, spamm_list_get_index(A_index.index, i));
  }

  for (i = 0; i < B_index.size; i++)
  {
    B_index.data[i] = spamm_hashtable_lookup(B_tier_hashtable, spamm_list_get_index(B_index.index, i));
  }

  spamm_timer_stop(timer);
  timer_string = spamm_timer_get_string(timer);
  printf("%s timer units\n", timer_string);
  free(timer_string);
#endif

#ifdef SPAMM_MULTIPLY_DOUBLE_CHECK
  for (i = 0, A_k_lookup_index = 0; i < A_index.size; i++)
  {
    if (spamm_list_get_norm(A_index.index, i) != A_index.data[i]->node_norm)
    {
      printf("[multiply] DOUBLECHECK norm mismatch in A_index[%u]\n", i);
      exit(1);
    }

    A_k = spamm_list_get_index(A_index.index, i) & MASK_2D_J;

    if (i == A_k_lookup.index[A_k_lookup_index])
    {
      first_A_k = spamm_list_get_index(A_index.index, i) & MASK_2D_J;
    }

    else if (i == A_k_lookup.index[A_k_lookup_index+1])
    {
      A_k_lookup_index++;
      first_A_k = spamm_list_get_index(A_index.index, i) & MASK_2D_J;
    }

    else
    {
      if (A_k != first_A_k)
      {
        printf("[multipy] DOUBLECHECK A_k_lookup incorrect\n");
        exit(1);
      }
    }

    if (i < A_index.size-1)
    {
      next_A_k = spamm_list_get_index(A_index.index, i+1) & MASK_2D_J;

      if (A_k == next_A_k) {
        if (A_index.data[i]->node_norm < A_index.data[i+1]->node_norm)
        {
          printf("[multiply] DOUBLECHECK norms in A_index are not sorted, norm[%u] = %e, norm[%u] = %e\n",
              i, A_index.data[i]->node_norm, i+1, A_index.data[i+1]->node_norm);
          exit(1);
        }
      }
    }
  }

  for (i = 0, B_k_lookup_index = 0; i < B_index.size; i++)
  {
    if (spamm_list_get_norm(B_index.index, i) != B_index.data[i]->node_norm)
    {
      printf("[multiply] DOUBLECHECK norm mismatch in B_index[%u]\n", i);
      exit(1);
    }

    B_k = spamm_list_get_index(B_index.index, i) & MASK_2D_I;

    if (i == B_k_lookup.index[B_k_lookup_index])
    {
      first_B_k = spamm_list_get_index(B_index.index, i) & MASK_2D_I;
    }

    else if (i == B_k_lookup.index[B_k_lookup_index+1])
    {
      B_k_lookup_index++;
      first_B_k = spamm_list_get_index(B_index.index, i) & MASK_2D_I;
    }

    else
    {
      if (B_k != first_B_k)
      {
        printf("[multipy] DOUBLECHECK B_k_lookup incorrect\n");
        exit(1);
      }
    }

    if (i < B_index.size-1)
    {
      next_B_k = spamm_list_get_index(B_index.index, i+1) & MASK_2D_I;

      if (B_k == next_B_k) {
        if (B_index.data[i]->node_norm < B_index.data[i+1]->node_norm)
        {
          printf("[multiply] DOUBLECHECK norms in B_index are not sorted, norm[%u] = %e, norm[%u] = %e\n",
              i, B_index.data[i]->node_norm, i+1, B_index.data[i+1]->node_norm);
          exit(1);
        }
      }
    }
  }
#endif

#ifdef SPAMM_MULTIPLY_CONVOLUTE
  /* Convolute by constructing product 3D index. */
  printf("[multiply] convolute... ");
  spamm_timer_start(timer);

  multiply_stream = spamm_allocate(sizeof(struct spamm_multiply_stream_t)
      *(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL));
  C_block_stream_index = spamm_allocate(sizeof(unsigned int)
      *(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL));

  /* Some tings to try:
   *
   * 1) Terminate the 2D index with a leading "1" bit, like Warren/Salmon to
   * indicate the width of the key and therefore its tier.
   */

#if SPAMM_MULTIPLY_CONVOLUTE_IMPLEMENTATION == 1
  stream_index = 0;
  number_dropped_blocks = 0;

  /* Loop over A. */
  for (A_k_lookup_index = 0, B_k_lookup_index = 0; A_k_lookup_index < A_k_lookup.size-1 && B_k_lookup_index < B_k_lookup.size-1; A_k_lookup_index++)
  {
    /* Get k value of A. */
    A_k = spamm_list_get_index(A_index.index, A_k_lookup.index[A_k_lookup_index]) & MASK_2D_J;

    /* Note that we don't increment i in the for() construct. */
    for (i = A_k_lookup.index[A_k_lookup_index]; i < A_k_lookup.index[A_k_lookup_index+1]; )
    {
      /* Get k value of B. */
      B_k = (spamm_list_get_index(B_index.index, B_k_lookup.index[B_k_lookup_index]) & MASK_2D_I) >> 1;

      /* Compare k values. */
      if (A_k > B_k)
      {
        /* Advance B in k. */
        B_k_lookup_index++;

        /* Possibly terminate. */
        if (B_k_lookup_index == B_k_lookup.size-1)
        {
          break;
        }
        continue;
      }

      else if (A_k < B_k)
      {
        continue;
      }

      /* Get reference to dense block of A. */
      A_data = A_index.data[i];

      /* Loop over subset of B with matching k. */
      for (j = B_k_lookup.index[B_k_lookup_index]; j < B_k_lookup.index[B_k_lookup_index+1]; j++)
      {
        /* Get reference to dense block of B. */
        B_data = B_index.data[j];

        /* Perform norm product and test whether to keep this term. */
        if (A_data->node_norm*B_data->node_norm <= tolerance)
        {
          break;
        }

        /* Get the linear 2D index of the C block. */
        convolution_index_2D = (spamm_list_get_index(A_index.index, i) & MASK_2D_I) |
          (spamm_list_get_index(B_index.index, j) & MASK_2D_J);

        /* Set references to matrix block in multiply stream. */
        multiply_stream[stream_index].A = A_data;
        multiply_stream[stream_index].B = B_data;

        /* Get reference to dense block of C. */
        C_data = spamm_hashtable_lookup(C_tier_hashtable, convolution_index_2D);

        if (C_data == NULL)
        {
          //printf("[%s:%i] creating new node of C\n", __FILE__, __LINE__);
          C_data = spamm_hashed_new_data(C->kernel_tier, convolution_index_2D, C->layout);
          spamm_hashtable_insert(C_tier_hashtable, convolution_index_2D, C_data);
        }

        multiply_stream[stream_index].C = C_data;

        /* Done with this stream element. */
        stream_index++;
      }

      /* Test how quickly we broke out in the previous loop over B. */
      if (j == B_k_lookup.index[B_k_lookup_index])
      {
        /* We never went past the first block in B. Since k segments are norm
         * sorted, we know that we can skip the rest of this k segment in A. */

        /* Increment dropped count. */
        number_dropped_blocks += (A_k_lookup.index[A_k_lookup_index+1]-i)*(B_k_lookup.index[B_k_lookup_index+1]-B_k_lookup.index[B_k_lookup_index]);

        /* Go to the next k-block in A. */
        A_k_lookup_index++;

        /* Possibly Terminate. */
        if (A_k_lookup_index == A_k_lookup.size-1)
        {
          break;
        }

        /* Set loop counter correctly. */
        i = A_k_lookup.index[A_k_lookup_index];

        /* Get k value of A. */
        A_k = spamm_list_get_index(A_index.index, A_k_lookup.index[A_k_lookup_index]) & MASK_2D_J;

        continue;
      }

      else
      {
        /* Increment dropped count. */
        number_dropped_blocks += B_k_lookup.index[B_k_lookup_index+1]-j;
      }

      /* Increment loop counter. */
      i++;
    }
  }
#elif SPAMM_MULTIPLY_CONVOLUTE_IMPLEMENTATION == 2
  unsigned int *A_index_array = spamm_list_get_index_address(A_index.index);
  unsigned int *B_index_array = spamm_list_get_index_address(B_index.index);

  unsigned int A_index_current[4];
  unsigned int B_index_current[4];
  unsigned int C_index_current[4];

  unsigned int A_index_k_current[4];
  unsigned int B_index_k_current[4];

  float A_norm_current[4];
  float B_norm_current[4];

  short int norm_product[4];

  short int early_termination = 0;

  stream_index = 0;
  number_dropped_blocks = 0;

  for (A_k_lookup_index = 0, B_k_lookup_index = 0; A_k_lookup_index < A_k_lookup.size && B_k_lookup_index < B_k_lookup.size; )
  {
    /* Match k-indices. */
    A_k = A_index_array[A_k_lookup.index[A_k_lookup_index]] & MASK_2D_J;
    B_k = (B_index_array[B_k_lookup.index[B_k_lookup_index]] & MASK_2D_I) >> 1;

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
      A_index_current[0] = A_index_array[i];
      A_index_current[1] = A_index_current[0];
      A_index_current[2] = A_index_current[0];
      A_index_current[3] = A_index_current[0];

      /* Load norm of A. */
      A_norm_current[0] = A_index.data[i]->node_norm;
      A_norm_current[1] = A_norm_current[0];
      A_norm_current[2] = A_norm_current[0];
      A_norm_current[3] = A_norm_current[0];

      early_termination = 0;
      for (j = B_k_lookup.index[B_k_lookup_index]; j < B_k_lookup.index[B_k_lookup_index+1]; j += 4)
      {
        /* Get index of B. */
        B_index_current[0] = B_index_array[j];
        B_index_current[1] = (j+1 < B_k_lookup.index[B_k_lookup_index+1] ? B_index_array[j+1] : 0);
        B_index_current[2] = (j+2 < B_k_lookup.index[B_k_lookup_index+1] ? B_index_array[j+2] : 0);
        B_index_current[3] = (j+3 < B_k_lookup.index[B_k_lookup_index+1] ? B_index_array[j+3] : 0);

        /* Load norm of B. */
        B_norm_current[0] = B_index.data[j]->node_norm;
        B_norm_current[1] = (j+1 < B_k_lookup.index[B_k_lookup_index+1] ? B_index.data[j+1]->node_norm : 0.0);
        B_norm_current[2] = (j+2 < B_k_lookup.index[B_k_lookup_index+1] ? B_index.data[j+2]->node_norm : 0.0);
        B_norm_current[3] = (j+3 < B_k_lookup.index[B_k_lookup_index+1] ? B_index.data[j+3]->node_norm : 0.0);

        /* Calculate norm products. */
        norm_product[0] = (A_norm_current[0]*B_norm_current[0] > tolerance) && (A_norm_current[0]*B_norm_current[0] > 0.0);
        norm_product[1] = (A_norm_current[1]*B_norm_current[1] > tolerance) && (A_norm_current[1]*B_norm_current[1] > 0.0);
        norm_product[2] = (A_norm_current[2]*B_norm_current[2] > tolerance) && (A_norm_current[2]*B_norm_current[2] > 0.0);
        norm_product[3] = (A_norm_current[3]*B_norm_current[3] > tolerance) && (A_norm_current[3]*B_norm_current[3] > 0.0);

        /* Calculate indices of C. */
        C_index_current[0] = (A_index_current[0] & MASK_2D_I) | (B_index_current[0] & MASK_2D_J);
        C_index_current[1] = (A_index_current[1] & MASK_2D_I) | (B_index_current[1] & MASK_2D_J);
        C_index_current[2] = (A_index_current[2] & MASK_2D_I) | (B_index_current[2] & MASK_2D_J);
        C_index_current[3] = (A_index_current[3] & MASK_2D_I) | (B_index_current[3] & MASK_2D_J);

        /* Append to stream or done. */
        if (norm_product[0])
        {
          multiply_stream[stream_index].A = A_index.data[i];
          multiply_stream[stream_index].B = B_index.data[j];
          multiply_stream[stream_index].C = spamm_hashtable_lookup(C_tier_hashtable, C_index_current[0]);
          stream_index++;
        }

        else
        {
          early_termination = 1;
          break;
        }

        if (norm_product[1])
        {
          multiply_stream[stream_index].A = A_index.data[i];
          multiply_stream[stream_index].B = B_index.data[j+1];
          multiply_stream[stream_index].C = spamm_hashtable_lookup(C_tier_hashtable, C_index_current[1]);
          stream_index++;
        }

        else { break; }

        if (norm_product[2])
        {
          multiply_stream[stream_index].A = A_index.data[i];
          multiply_stream[stream_index].B = B_index.data[j+2];
          multiply_stream[stream_index].C = spamm_hashtable_lookup(C_tier_hashtable, C_index_current[2]);
          stream_index++;
        }

        else { break; }

        if (norm_product[3])
        {
          multiply_stream[stream_index].A = A_index.data[i];
          multiply_stream[stream_index].B = B_index.data[j+3];
          multiply_stream[stream_index].C = spamm_hashtable_lookup(C_tier_hashtable, C_index_current[3]);
          stream_index++;
        }

        else { break; }
      }

      /* Check whether we had at least one block in B that made it passed the
       * SpAMM condition. */
      if (early_termination == 1)
      {
        break;
      }
    }

    A_k_lookup_index++;
  }
#elif SPAMM_MULTIPLY_CONVOLUTE_IMPLEMENTATION == 3
  unsigned int *A_index_array = spamm_list_get_index_address(A_index.index);
  unsigned int *B_index_array = spamm_list_get_index_address(B_index.index);

  unsigned int *multiply_stream_C_index = spamm_allocate(sizeof(unsigned int)
      *(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL));

  unsigned int A_index_current[4];
  unsigned int B_index_current[4];
  unsigned int C_index_current[4];

  unsigned int A_index_k_current[4];
  unsigned int B_index_k_current[4];

  float A_norm_current[4];
  float B_norm_current[4];

  short int norm_product[4];

  short int early_termination = 0;

  unsigned int last_C_index;
  struct spamm_hashed_data_t *last_C;

  stream_index = 0;
  number_dropped_blocks = 0;

  for (A_k_lookup_index = 0, B_k_lookup_index = 0; A_k_lookup_index < A_k_lookup.size && B_k_lookup_index < B_k_lookup.size; )
  {
    /* Match k-indices. */
    A_k = A_index_array[A_k_lookup.index[A_k_lookup_index]] & MASK_2D_J;
    B_k = (B_index_array[B_k_lookup.index[B_k_lookup_index]] & MASK_2D_I) >> 1;

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
      A_index_current[0] = A_index_array[i];
      A_index_current[1] = A_index_current[0];
      A_index_current[2] = A_index_current[0];
      A_index_current[3] = A_index_current[0];

      /* Load norm of A. */
      A_norm_current[0] = A_index.data[i]->node_norm;
      A_norm_current[1] = A_norm_current[0];
      A_norm_current[2] = A_norm_current[0];
      A_norm_current[3] = A_norm_current[0];

      early_termination = 0;
      for (j = B_k_lookup.index[B_k_lookup_index]; j < B_k_lookup.index[B_k_lookup_index+1]; j += 4)
      {
        /* Get index of B. */
        B_index_current[0] = B_index_array[j];
        B_index_current[1] = (j+1 < B_k_lookup.index[B_k_lookup_index+1] ? B_index_array[j+1] : 0);
        B_index_current[2] = (j+2 < B_k_lookup.index[B_k_lookup_index+1] ? B_index_array[j+2] : 0);
        B_index_current[3] = (j+3 < B_k_lookup.index[B_k_lookup_index+1] ? B_index_array[j+3] : 0);

        /* Load norm of B. */
        B_norm_current[0] = B_index.data[j]->node_norm;
        B_norm_current[1] = (j+1 < B_k_lookup.index[B_k_lookup_index+1] ? B_index.data[j+1]->node_norm : 0.0);
        B_norm_current[2] = (j+2 < B_k_lookup.index[B_k_lookup_index+1] ? B_index.data[j+2]->node_norm : 0.0);
        B_norm_current[3] = (j+3 < B_k_lookup.index[B_k_lookup_index+1] ? B_index.data[j+3]->node_norm : 0.0);

        /* Calculate norm products. */
        norm_product[0] = (A_norm_current[0]*B_norm_current[0] > tolerance) && (A_norm_current[0]*B_norm_current[0] > 0.0);
        norm_product[1] = (A_norm_current[1]*B_norm_current[1] > tolerance) && (A_norm_current[1]*B_norm_current[1] > 0.0);
        norm_product[2] = (A_norm_current[2]*B_norm_current[2] > tolerance) && (A_norm_current[2]*B_norm_current[2] > 0.0);
        norm_product[3] = (A_norm_current[3]*B_norm_current[3] > tolerance) && (A_norm_current[3]*B_norm_current[3] > 0.0);

        /* Calculate indices of C. */
        C_index_current[0] = (A_index_current[0] & MASK_2D_I) | (B_index_current[0] & MASK_2D_J);
        C_index_current[1] = (A_index_current[1] & MASK_2D_I) | (B_index_current[1] & MASK_2D_J);
        C_index_current[2] = (A_index_current[2] & MASK_2D_I) | (B_index_current[2] & MASK_2D_J);
        C_index_current[3] = (A_index_current[3] & MASK_2D_I) | (B_index_current[3] & MASK_2D_J);

        /* Append to stream or done. */
        if (norm_product[0])
        {
          multiply_stream[stream_index].A = A_index.data[i];
          multiply_stream[stream_index].B = B_index.data[j];
          multiply_stream_C_index[stream_index] = C_index_current[0];
          stream_index++;
        }

        else
        {
          early_termination = 1;
          break;
        }

        if (norm_product[1])
        {
          multiply_stream[stream_index].A = A_index.data[i];
          multiply_stream[stream_index].B = B_index.data[j+1];
          multiply_stream_C_index[stream_index] = C_index_current[1];
          stream_index++;
        }

        else { break; }

        if (norm_product[2])
        {
          multiply_stream[stream_index].A = A_index.data[i];
          multiply_stream[stream_index].B = B_index.data[j+2];
          multiply_stream_C_index[stream_index] = C_index_current[2];
          stream_index++;
        }

        else { break; }

        if (norm_product[3])
        {
          multiply_stream[stream_index].A = A_index.data[i];
          multiply_stream[stream_index].B = B_index.data[j+3];
          multiply_stream_C_index[stream_index] = C_index_current[3];
          stream_index++;
        }

        else { break; }
      }

      /* Check whether we had at least one block in B that made it passed the
       * SpAMM condition. */
      if (early_termination == 1)
      {
        break;
      }
    }

    A_k_lookup_index++;
  }

  /* Sort multiply_stream and lookup C blocks. */
  spamm_multiply_C_index_sort(multiply_stream_C_index, multiply_stream, stream_index);

  last_C_index = multiply_stream_C_index[0];
  last_C = spamm_hashtable_lookup(C_tier_hashtable, last_C_index);
  for (i = 0; i < stream_index; i++)
  {
    if (multiply_stream_C_index[i] != last_C_index)
    {
      last_C_index = multiply_stream_C_index[i];
      last_C = spamm_hashtable_lookup(C_tier_hashtable, last_C_index);
    }

    /* Assign C block. */
    multiply_stream[i].C = last_C;
  }
#elif SPAMM_MULTIPLY_CONVOLUTE_IMPLEMENTATION == 4
  unsigned int *A_index_array = spamm_list_get_index_address(A_index.index);
  unsigned int *B_index_array = spamm_list_get_index_address(B_index.index);

  unsigned int *multiply_stream_C_index = spamm_allocate(sizeof(unsigned int)
      *(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL));

  unsigned int A_index_current[4];
  unsigned int B_index_current[4];
  unsigned int C_index_current[4];

  float A_norm_current[4];
  float B_norm_current[4];

  short int norm_product[4];

  short int early_termination = 0;

  unsigned int last_C_index;
  struct spamm_hashed_data_t *last_C;

  stream_index = 0;
  number_dropped_blocks = 0;

  for (A_k_lookup_index = 0, B_k_lookup_index = 0; A_k_lookup_index+1 < A_k_lookup.size && B_k_lookup_index+1 < B_k_lookup.size; )
  {
    /* Match k-indices. */
    A_k = A_index_array[A_k_lookup.index[A_k_lookup_index]] & MASK_2D_J;
    B_k = (B_index_array[B_k_lookup.index[B_k_lookup_index]] & MASK_2D_I) >> 1;

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
      A_index_current[0] = A_index_array[i];
      A_index_current[1] = A_index_current[0];
      A_index_current[2] = A_index_current[0];
      A_index_current[3] = A_index_current[0];

      /* Load norm of A. */
      A_norm_current[0] = A_index.data[i]->node_norm;
      A_norm_current[1] = A_norm_current[0];
      A_norm_current[2] = A_norm_current[0];
      A_norm_current[3] = A_norm_current[0];

      for (j = B_k_lookup.index[B_k_lookup_index], early_termination = 1; j < B_k_lookup.index[B_k_lookup_index+1]; j += 4)
      {
        /* Get index of B. */
        B_index_current[0] = B_index_array[j];
        B_index_current[1] = (j+1 < B_k_lookup.index[B_k_lookup_index+1] ? B_index_array[j+1] : 0);
        B_index_current[2] = (j+2 < B_k_lookup.index[B_k_lookup_index+1] ? B_index_array[j+2] : 0);
        B_index_current[3] = (j+3 < B_k_lookup.index[B_k_lookup_index+1] ? B_index_array[j+3] : 0);

        /* Load norm of B. */
        B_norm_current[0] = B_index.data[j]->node_norm;
        B_norm_current[1] = (j+1 < B_k_lookup.index[B_k_lookup_index+1] ? B_index.data[j+1]->node_norm : 0.0);
        B_norm_current[2] = (j+2 < B_k_lookup.index[B_k_lookup_index+1] ? B_index.data[j+2]->node_norm : 0.0);
        B_norm_current[3] = (j+3 < B_k_lookup.index[B_k_lookup_index+1] ? B_index.data[j+3]->node_norm : 0.0);

        /* Calculate norm products. */
        norm_product[0] = (A_norm_current[0]*B_norm_current[0] > tolerance);
        norm_product[1] = (A_norm_current[1]*B_norm_current[1] > tolerance);
        norm_product[2] = (A_norm_current[2]*B_norm_current[2] > tolerance);
        norm_product[3] = (A_norm_current[3]*B_norm_current[3] > tolerance);

        /* Calculate indices of C. */
        C_index_current[0] = (A_index_current[0] & MASK_2D_I) | (B_index_current[0] & MASK_2D_J);
        C_index_current[1] = (A_index_current[1] & MASK_2D_I) | (B_index_current[1] & MASK_2D_J);
        C_index_current[2] = (A_index_current[2] & MASK_2D_I) | (B_index_current[2] & MASK_2D_J);
        C_index_current[3] = (A_index_current[3] & MASK_2D_I) | (B_index_current[3] & MASK_2D_J);

        /* Append to stream or done. */
        if (norm_product[0])
        {
          multiply_stream[stream_index].A = A_index.data[i];
          multiply_stream[stream_index].B = B_index.data[j];
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

        if (norm_product[1])
        {
          multiply_stream[stream_index].A = A_index.data[i];
          multiply_stream[stream_index].B = B_index.data[j+1];
          multiply_stream_C_index[stream_index] = C_index_current[1];
          stream_index++;
        }

        else
        {
          number_dropped_blocks += B_k_lookup.index[B_k_lookup_index+1]-(j+1);
          break;
        }

        if (norm_product[2])
        {
          multiply_stream[stream_index].A = A_index.data[i];
          multiply_stream[stream_index].B = B_index.data[j+2];
          multiply_stream_C_index[stream_index] = C_index_current[2];
          stream_index++;
        }

        else
        {
          number_dropped_blocks += B_k_lookup.index[B_k_lookup_index+1]-(j+2);
          break;
        }

        if (norm_product[3])
        {
          multiply_stream[stream_index].A = A_index.data[i];
          multiply_stream[stream_index].B = B_index.data[j+3];
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

    /* Assign C block. */
    multiply_stream[i].C = last_C;
  }
#elif SPAMM_MULTIPLY_CONVOLUTE_IMPLEMENTATION == 5
  unsigned int *A_index_array = spamm_list_get_index_address(A_index.index);
  unsigned int *B_index_array_original = spamm_list_get_index_address(B_index.index);

  unsigned int *multiply_stream_C_index = spamm_allocate(sizeof(unsigned int)
      *(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL));

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

  A_norm_array = spamm_allocate(A_index.size*sizeof(float));
  for (i = 0; i < A_index.size; i++)
  {
    A_norm_array[i] = A_index.data[i]->node_norm;
  }

  /* Copy index and norm arrays. This needs to be done to create properly
   * aligned arrays. */
  B_k_lookup_array = spamm_allocate(B_k_lookup.size*sizeof(unsigned int));
  B_norm_array = spamm_allocate(B_index.size*sizeof(float)+B_k_lookup.size*16);
  B_index_array = spamm_allocate(B_index.size*sizeof(unsigned int)+B_k_lookup.size*16);
  B_index_translation = spamm_allocate(B_index.size*sizeof(unsigned int)+B_k_lookup.size*16);

  for (B_k_lookup_index = 0, j = 0; B_k_lookup_index+1 < B_k_lookup.size; B_k_lookup_index++)
  {
    B_k_lookup_array[B_k_lookup_index] = j;

    /* Copy k-block data. */
    for (i = B_k_lookup.index[B_k_lookup_index]; i < B_k_lookup.index[B_k_lookup_index+1]; i++)
    {
      B_norm_array[j] = B_index.data[i]->node_norm;
      B_index_array[j] = B_index_array_original[i];
      B_index_translation[j] = i;
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

#endif

  /* Check. */
  if (stream_index > (A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL))
  {
    printf("multiply stream has too many elements, has %u but is only dimensioned for %u\n", stream_index,
        (A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL));
    exit(1);
  }

#ifdef SPAMM_MULTIPLY_DOUBLE_CHECK
  /* Print multiply stream for debugging. */
  printf("[multiply] multiply stream = {\n");
  for (i = 0; i < stream_index; i++)
  {
    printf("  [ %u(%e) %u(%e) -> %u(%e) ] : %u\n",
        multiply_stream[i].A->index_2D, multiply_stream[i].A->node_norm,
        multiply_stream[i].B->index_2D, multiply_stream[i].B->node_norm,
        multiply_stream[i].C->index_2D, multiply_stream[i].A->node_norm*multiply_stream[i].B->node_norm, i);
  }
  printf("}\n");
#endif

  spamm_timer_stop(timer);
  timer_string = spamm_timer_get_string(timer);
  printf("%s timer units\n", timer_string);
  free(timer_string);
  printf("[multiply] dropped %u blocks, placed %u blocks into stream, total of %u blocks\n", number_dropped_blocks, stream_index, number_dropped_blocks+stream_index);
#endif

#ifdef SPAMM_MULTIPLY_FREE
  /* Free memory. */
  printf("[multiply] free memory... ");
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
  printf("%s timer units\n", timer_string);
  free(timer_string);
#endif

#ifdef SPAMM_MULTIPLY_STREAM
  /* Call stream product. */
  printf("[multiply] stream multiply ");
  spamm_timer_start(timer);

  printf("(%s)... ", spamm_kernel_get_name(kernel));
  switch (kernel)
  {
    case kernel_external_sgemm:
      spamm_stream_external_sgemm(stream_index, alpha, tolerance, multiply_stream, 1);
      break;

    case kernel_stream_NULL:
      spamm_stream_external_sgemm(stream_index, alpha, tolerance, multiply_stream, 0);
      break;

    case kernel_standard_SSE:
      spamm_stream_kernel_SSE(stream_index, alpha, tolerance, multiply_stream);
      break;

    case kernel_standard_SSE4_1:
      spamm_stream_kernel_SSE4_1(stream_index, alpha, tolerance, multiply_stream);
      break;

    default:
      printf("[multiply] unknown kernel... ");
      exit(1);
      break;
  }

  spamm_timer_stop(timer);
  timer_string = spamm_timer_get_string(timer);
  printf("%s timer units\n", timer_string);
  free(timer_string);
#endif

#ifdef SPAMM_MULTIPLY_UPDATE_NORM
  printf("[multiply] updating C tree norms... ");
  spamm_timer_start(timer);

  spamm_hashed_norm_update(C);
  spamm_construct_tree(C);

  spamm_timer_stop(timer);
  timer_string = spamm_timer_get_string(timer);
  printf("%s timer units\n", timer_string);
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
 */
void
spamm_recursive_multiply_matrix (const float tolerance,
    const float alpha,
    struct spamm_recursive_node_t *node_A,
    struct spamm_recursive_node_t *node_B,
    struct spamm_recursive_node_t **node_C,
    struct spamm_timer_t *timer,
    void (*sgemm) (char *, char *, int *, int *, int *, float *, float *, int *, float *, int *, float *, float *, int *),
    unsigned int *number_products)
{
  float beta = 1.0;
  int i, j, k;

  if (node_A == NULL || node_B == NULL) { return; }

  /* We have to allocate a new C block a tier up. */
  if (*node_C == NULL)
  {
    printf("[%s:%i] node_C should not be NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  /* Multiply matrix blocks. */
  if (node_A->data != NULL && node_B->data != NULL)
  {
    if (node_A->norm*node_B->norm > tolerance)
    {
      if ((*node_C)->data == NULL)
      {
        (*node_C)->data = calloc((*node_C)->N_contiguous*(*node_C)->N_contiguous, sizeof(float));
      }
      if (sgemm != NULL)
      {
        sgemm(
            "N", /* TRANSA */
            "N", /* TRANSB */
            (int*) &(node_A->N_contiguous), /* M */
            (int*) &(node_A->N_contiguous), /* N */
            (int*) &(node_A->N_contiguous), /* K */
            (float*) &alpha, /* alpha */
            node_A->data, /* A */
            (int*) &node_A->N_contiguous, /* LDA */
            node_B->data, /* B */
            (int*) &node_A->N_contiguous, /* LDB */
            (float*) &beta, /* beta */
            (*node_C)->data, /* C */
            (int*) &node_A->N_contiguous /* LDC */
            );
        (*number_products)++;
      }
    }
  }

  else
  {
    /* Recurse. */
    for (i = 0; i < 2; i++) {
      for (j = 0; j < 2; j++) {
        for (k = 0; k < 2; k++)
        {
          if (node_A->child[spamm_index_row_major(i, k, 2, 2)] != NULL && node_B->child[spamm_index_row_major(k, j, 2, 2)] != NULL)
          {
            if (node_A->child[spamm_index_row_major(i, k, 2, 2)]->norm *
                node_B->child[spamm_index_row_major(k, j, 2, 2)]->norm > tolerance)
            {
              /* Create a new C node if necessary. */
              if ((*node_C)->child[spamm_index_row_major(i, j, 2, 2)] == NULL)
              {
                (*node_C)->child[spamm_index_row_major(i, j, 2, 2)] = spamm_recursive_new_node((*node_C)->tier+1, (*node_C)->N_contiguous);
                (*node_C)->child[spamm_index_row_major(i, j, 2, 2)]->M_lower = (*node_C)->M_lower+((*node_C)->M_upper-(*node_C)->M_lower)/2*i;
                (*node_C)->child[spamm_index_row_major(i, j, 2, 2)]->M_upper = (*node_C)->M_lower+((*node_C)->M_upper-(*node_C)->M_lower)/2*(i+1);
                (*node_C)->child[spamm_index_row_major(i, j, 2, 2)]->N_lower = (*node_C)->N_lower+((*node_C)->N_upper-(*node_C)->N_lower)/2*j;
                (*node_C)->child[spamm_index_row_major(i, j, 2, 2)]->N_upper = (*node_C)->N_lower+((*node_C)->N_upper-(*node_C)->N_lower)/2*(j+1);
              }

              spamm_recursive_multiply_matrix(tolerance,
                  alpha,
                  node_A->child[spamm_index_row_major(i, k, 2, 2)],
                  node_B->child[spamm_index_row_major(k, j, 2, 2)],
                  &((*node_C)->child[spamm_index_row_major(i, j, 2, 2)]),
                  timer,
                  sgemm,
                  number_products);
            }
          }
        }
      }
    }
  }
}

/** @brief Multiply a node by a scalar.
 *
 * The node and its children nodes are multiplied by the scalar.
 *
 * @param beta The scalar.
 * @param node The node.
 */
void
spamm_recursive_multiply_scalar (const float beta, struct spamm_recursive_node_t *node)
{
  int i;

  if (node == NULL) { return; }

  if (node->data != NULL)
  {
    for (i = 0; i < node->N_contiguous*node->N_contiguous; i++)
    {
      node->data[i] *= beta;
    }
  }

  else
  {
    for (i = 0; i < 4; i++)
    {
      spamm_recursive_multiply_scalar(beta, node->child[i]);
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
 * @param sgemm The external sgemm function to use.
 */
void
spamm_recursive_multiply (const float tolerance,
    const float alpha,
    struct spamm_recursive_t *A,
    struct spamm_recursive_t *B,
    const float beta,
    struct spamm_recursive_t *C,
    struct spamm_timer_t *timer,
    void (*sgemm) (char *, char *, int *, int *, int *, float *, float *, int *, float *, int *, float *, float *, int *),
    unsigned int *number_products)
{
  if (A == NULL)
  {
    printf("[%s:%i] A is NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  if (B == NULL)
  {
    printf("[%s:%i] B is NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  if (C == NULL)
  {
    printf("[%s:%i] C is NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  if (A->N_contiguous != B->N_contiguous)
  {
    printf("[%s:%i] A->N_contiguous != B->N_contiguous\n", __FILE__, __LINE__);
    exit(1);
  }

  if (A->N_contiguous != C->N_contiguous)
  {
    printf("[%s:%i] A->N_contiguous != C->N_contiguous\n", __FILE__, __LINE__);
    exit(1);
  }

  /* Multiply C by beta. */
  spamm_recursive_multiply_scalar(beta, C->root);

  /* Multiply A and B. */
  if (A->root != NULL && B->root != NULL && C->root == NULL)
  {
    C->root = spamm_recursive_new_node(0, C->N_contiguous);
    C->root->M_lower = 0;
    C->root->M_upper = A->root->M_upper;
    C->root->N_lower = 0;
    C->root->N_upper = A->root->N_upper;
  }

  spamm_recursive_multiply_matrix(tolerance, alpha, A->root, B->root, &(C->root), timer, sgemm, number_products);
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
    void (*sgemm) (char *, char *, int *, int *, int *, float *, float *, int *, float *, int *, float *, float *, int *),
    const enum spamm_kernel_t kernel,
    unsigned int *number_products)
{
  if (A == NULL)
  {
    printf("[%s:%i] A is NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  if (B == NULL)
  {
    printf("[%s:%i] B is NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  if (C == NULL)
  {
    printf("[%s:%i] C is NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  if (A->N_contiguous != B->N_contiguous)
  {
    printf("[%s:%i] A->N_contiguous != B->N_contiguous\n", __FILE__, __LINE__);
    exit(1);
  }

  if (A->N_contiguous != C->N_contiguous)
  {
    printf("[%s:%i] A->N_contiguous != C->N_contiguous\n", __FILE__, __LINE__);
    exit(1);
  }

  /* Multiply C by beta. */
  //spamm_recursive_multiply_scalar(beta, C->root);

  /* Multiply A and B. */
  //if (A->root != NULL && B->root != NULL && C->root == NULL)
  //{
  //  C->root = spamm_recursive_new_node(0, C->N_contiguous);
  //  C->root->M_lower = 0;
  //  C->root->M_upper = A->root->M_upper;
  //  C->root->N_lower = 0;
  //  C->root->N_upper = A->root->N_upper;
  //}

  //spamm_recursive_multiply_matrix(tolerance, alpha, A->root, B->root, &(C->root), timer, sgemm, number_products);
}

/** Multiply three matrices, i.e. \f$ D = \alpha A \times B \times C + \beta D\f$.
 *
 * @param tolerance The SpAMM tolerance of this product.
 * @param alpha The paramater \f$\alpha\f$.
 * @param A The matrix \f$A\f$.
 * @param B The matrix \f$B\f$.
 * @param C The matrix \f$C\f$.
 * @param beta The paramater \f$\beta\f$.
 * @param D The matrix \f$D\f$.
 * @param timer The timer to use.
 * @param sgemm The external sgemm function to use.
 */
void
spamm_recursive_multiply_3_matrix (const float tolerance,
    const float alpha,
    struct spamm_recursive_node_t *node_A,
    struct spamm_recursive_node_t *node_B,
    struct spamm_recursive_node_t *node_C,
    struct spamm_recursive_node_t **node_D,
    struct spamm_timer_t *timer,
    void (*sgemm) (char *, char *, int *, int *, int *, float *, float *, int *, float *, int *, float *, float *, int *),
    unsigned int *number_products)
{
  float beta = 1.0;
  int i, j, k;

  if (node_A == NULL || node_B == NULL || node_C == NULL) { return; }

  /* We have to allocate a new D block a tier up. */
  if (*node_D == NULL)
  {
    printf("[%s:%i] node_D should not be NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  /* Multiply matrix blocks. */
  if (node_A->data != NULL && node_B->data != NULL && node_C->data != NULL)
  {
    if (node_A->norm*node_B->norm*node_C->norm > tolerance)
    {
      if ((*node_D)->data == NULL)
      {
        (*node_D)->data = calloc((*node_D)->N_contiguous*(*node_D)->N_contiguous, sizeof(float));
      }
      if (sgemm != NULL)
      {
        sgemm("N", "N", (int*) &(node_A->N_contiguous), (int*) &(node_A->N_contiguous), (int*) &(node_A->N_contiguous),
            (float*) &alpha, node_A->data, (int*) &node_A->N_contiguous, node_B->data,
            (int*) &node_A->N_contiguous, (float*) &beta, (*node_D)->data, (int*) &node_A->N_contiguous);
        (*number_products)++;
      }
    }
  }

//  else
//  {
//    /* Recurse. */
//    for (i = 0; i < 2; i++) {
//      for (j = 0; j < 2; j++) {
//        for (k = 0; k < 2; k++)
//        {
//          if (node_A->child[spamm_index_row_major(i, k, 2, 2)] != NULL && node_B->child[spamm_index_row_major(k, j, 2, 2)] != NULL)
//          {
//            if (node_A->child[spamm_index_row_major(i, k, 2, 2)]->norm *
//                node_B->child[spamm_index_row_major(k, j, 2, 2)]->norm > tolerance)
//            {
//              /* Create a new C node if necessary. */
//              if ((*node_C)->child[spamm_index_row_major(i, j, 2, 2)] == NULL)
//              {
//                (*node_C)->child[spamm_index_row_major(i, j, 2, 2)] = spamm_recursive_new_node((*node_C)->tier+1, (*node_C)->blocksize);
//                (*node_C)->child[spamm_index_row_major(i, j, 2, 2)]->M_lower = (*node_C)->M_lower+((*node_C)->M_upper-(*node_C)->M_lower)/2*i;
//                (*node_C)->child[spamm_index_row_major(i, j, 2, 2)]->M_upper = (*node_C)->M_lower+((*node_C)->M_upper-(*node_C)->M_lower)/2*(i+1);
//                (*node_C)->child[spamm_index_row_major(i, j, 2, 2)]->N_lower = (*node_C)->N_lower+((*node_C)->N_upper-(*node_C)->N_lower)/2*j;
//                (*node_C)->child[spamm_index_row_major(i, j, 2, 2)]->N_upper = (*node_C)->N_lower+((*node_C)->N_upper-(*node_C)->N_lower)/2*(j+1);
//              }
//
//              spamm_recursive_multiply_matrix(tolerance,
//                  alpha,
//                  node_A->child[spamm_index_row_major(i, k, 2, 2)],
//                  node_B->child[spamm_index_row_major(k, j, 2, 2)],
//                  &((*node_C)->child[spamm_index_row_major(i, j, 2, 2)]),
//                  timer,
//                  sgemm,
//                  number_products);
//            }
//          }
//        }
//      }
//    }
//  }
}

/** Multiply three matrices, i.e. \f$ D = \alpha A \times B \times C + \beta D\f$.
 *
 * @param tolerance The SpAMM tolerance of this product.
 * @param alpha The paramater \f$\alpha\f$.
 * @param A The matrix \f$A\f$.
 * @param B The matrix \f$B\f$.
 * @param B The matrix \f$C\f$.
 * @param beta The paramater \f$\beta\f$.
 * @param C The matrix \f$D\f$.
 * @param timer The timer to use.
 * @param sgemm The external sgemm function to use.
 */
void
spamm_recursive_multiply_3 (const float tolerance,
    const float alpha,
    struct spamm_recursive_t *A,
    struct spamm_recursive_t *B,
    struct spamm_recursive_t *C,
    const float beta,
    struct spamm_recursive_t *D,
    struct spamm_timer_t *timer,
    void (*sgemm) (),
    unsigned int *number_products)
{
  if (A == NULL)
  {
    printf("[%s:%i] A is NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  if (B == NULL)
  {
    printf("[%s:%i] B is NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  if (C == NULL)
  {
    printf("[%s:%i] C is NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  if (D == NULL)
  {
    printf("[%s:%i] D is NULL\n", __FILE__, __LINE__);
    exit(1);
  }

  if (A->N_contiguous != B->N_contiguous)
  {
    printf("[%s:%i] A->N_contiguous != B->N_contiguous\n", __FILE__, __LINE__);
    exit(1);
  }

  if (A->N_contiguous != C->N_contiguous)
  {
    printf("[%s:%i] A->N_contiguous != C->N_contiguous\n", __FILE__, __LINE__);
    exit(1);
  }

  if (A->N_contiguous != D->N_contiguous)
  {
    printf("[%s:%i] A->N_contiguous != D->N_contiguous\n", __FILE__, __LINE__);
    exit(1);
  }

  /* Multiply C by beta. */
  spamm_recursive_multiply_scalar(beta, D->root);

  /* Multiply A and B. */
  if (A->root != NULL && B->root != NULL && C->root == NULL && D->root != NULL)
  {
    D->root = spamm_recursive_new_node(0, D->N_contiguous);
    D->root->M_lower = 0;
    D->root->M_upper = A->root->M_upper;
    D->root->N_lower = 0;
    D->root->N_upper = A->root->N_upper;
  }

  spamm_recursive_multiply_3_matrix(tolerance, alpha, A->root, B->root, C->root, &(C->root), timer, sgemm, number_products);
}
