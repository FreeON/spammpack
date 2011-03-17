#include "spamm_config.h"
#include "spamm.h"
#include <assert.h>
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
  struct spamm_data_t **data;
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
 * @param value A pointer to the spamm_data_t matrix node.
 * @param user_data The scalar \f$\beta\f$ that multiplies the matrix.
 */
void
spamm_multiply_beta_block (unsigned int index, void *value, void *user_data)
{
  struct spamm_data_t *data = value;
  float *beta = user_data;
  unsigned int i;
  short j;

#ifdef HAVE_SSE_DISABLED
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
spamm_multiply_beta (const float beta, struct spamm_t *A)
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
  struct spamm_data_t *temp_node;
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

/** Multiply to matrices, i.e. \f$ C = \alpha A \times B + \beta C\f$.
 *
 * @param tolerance The SpAMM tolerance of this product.
 * @param alpha The paramater \f$\alpha\f$.
 * @param A The matrix \f$A\f$.
 * @param B The matrix \f$B\f$.
 * @param beta The paramater \f$\beta\f$.
 * @param C The matrix \f$C\f$.
 * @param timer_type The timer to use.
 */
void
spamm_multiply (const float tolerance,
    const float alpha, struct spamm_t *A, struct spamm_t *B,
    const float beta, struct spamm_t *C,
    const enum spamm_timer_type_t timer_type)
{
  struct spamm_hashtable_t *A_tier_hashtable;
  struct spamm_hashtable_t *B_tier_hashtable;
  struct spamm_hashtable_t *C_tier_hashtable;

  struct spamm_multiply_index_list_t A_index;
  struct spamm_multiply_index_list_t B_index;

  struct spamm_data_t *A_data;
  struct spamm_data_t *B_data;
  struct spamm_data_t *C_data;

  unsigned int i, j, k, k_check;
  unsigned int index;
  unsigned int convolution_index_2D;
  unsigned int A_k_lookup_index;
  unsigned int B_k_lookup_index;
  unsigned int A_k, B_k;

  struct spamm_multiply_k_lookup_t A_k_lookup;
  struct spamm_multiply_k_lookup_t B_k_lookup;

  struct spamm_multiply_stream_t *multiply_stream;
  unsigned int stream_index;
  unsigned int number_dropped_blocks;

  unsigned int *C_block_stream_index;

  unsigned int number_products = 0;

  struct spamm_timer_t *timer = spamm_timer_new(timer_type);

  char timer_info_string[2000];

  unsigned long long beta_timer;
  unsigned long long sort_timer;
  unsigned long long k_lookuptable_timer;
  unsigned long long copy_timer;
  unsigned long long copy_3D_timer;
  unsigned long long free_timer;
  unsigned long long expand_timer;
  unsigned long long convolute_timer;
  unsigned long long stream_timer;

  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);

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
  beta_timer = spamm_timer_get(timer);
  printf("%llu timer units\n", beta_timer);
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

  spamm_list_sort_index(A_index.index, spamm_list_compare_index_column);
  spamm_list_sort_index(B_index.index, spamm_list_compare_index_row);

  printf("len(A) = %u, len(B) = %u, ", spamm_list_length(A_index.index), spamm_list_length(B_index.index));

  spamm_timer_stop(timer);
  sort_timer = spamm_timer_get(timer);
  printf("%llu timer units\n", sort_timer);
#endif

#ifdef SPAMM_MULTIPLY_K_LOOKUP
  /* Create a lookup table for the start of a particular k index in the sorted
   * arrays. */
  printf("[multiply] creating k lookup tables... ");
  spamm_timer_start(timer);

  A_k_lookup.index = (unsigned int*) malloc(sizeof(unsigned int)*(A->N_padded/SPAMM_N_KERNEL+1));
  B_k_lookup.index = (unsigned int*) malloc(sizeof(unsigned int)*(B->N_padded/SPAMM_N_KERNEL+1));

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
  k_lookuptable_timer = spamm_timer_get(timer);
  printf("%llu timer units\n", k_lookuptable_timer);
#endif

#ifdef SPAMM_MULTIPLY_SORT_NORM
  /* Sort 2D indices on norms. */
  printf("[multiply] sorting A and B on norms... ");
  spamm_timer_start(timer);

  for (i = 0; i < A_k_lookup.size-1; i++)
  {
    spamm_list_sort_norm(A_index.index, A_k_lookup.index[i], A_k_lookup.index[i+1]);
  }

  //spamm_list_print(A_index.index);

  for (i = 0; i < B_k_lookup.size-1; i++)
  {
    spamm_list_sort_norm(B_index.index, B_k_lookup.index[i], B_k_lookup.index[i+1]);
  }

  spamm_timer_stop(timer);
  sort_timer = spamm_timer_get(timer);
  printf("%llu timer units\n", sort_timer);
#endif

#ifdef SPAMM_MULTIPLY_COPY_INDICES
  /* Copy sorted indices to array for quick access. */
  printf("[multiply] copying indices to array... ");
  spamm_timer_start(timer);

  A_index.size = spamm_list_length(A_index.index);
  A_index.data = (struct spamm_data_t**) malloc(sizeof(struct spamm_data_t*)*spamm_list_length(A_index.index));

  B_index.size = spamm_list_length(B_index.index);
  B_index.data = (struct spamm_data_t**) malloc(sizeof(struct spamm_data_t*)*spamm_list_length(B_index.index));

  printf("len(A_index) = %u, len(B_index) = %u, ", A_index.size, B_index.size);

  spamm_timer_stop(timer);
  copy_timer = spamm_timer_get(timer);
  printf("%llu timer units\n", copy_timer);
#endif

#ifdef SPAMM_MULTIPLY_COPY_3D
  /* Copy appropriate 3D convolution index to arrays. */
  printf("[multiply] copying 3D convolution index to arrays and referencing dense blocks... ");
  spamm_timer_start(timer);

  for (i = 0; i < A_index.size; i++)
  {
    A_index.data[i] = spamm_hashtable_lookup(A_tier_hashtable, spamm_list_get_index(A_index.index, i));
  }

  for (i = 0; i < B_index.size; i++)
  {
    B_index.data[i] = spamm_hashtable_lookup(B_tier_hashtable, spamm_list_get_index(B_index.index, i));
  }

  spamm_timer_stop(timer);
  copy_3D_timer = spamm_timer_get(timer);
  printf("%llu timer units\n", copy_3D_timer);
#endif

#ifdef SPAMM_MULTIPLY_CONVOLUTE
  /* Convolute by constructing product 3D index. */
  printf("[multiply] convolute... ");
  spamm_timer_start(timer);

  stream_index = 0;
  number_dropped_blocks = 0;
  multiply_stream = (struct spamm_multiply_stream_t*) malloc(sizeof(struct spamm_multiply_stream_t)
      *(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL));
  C_block_stream_index = (unsigned int*) malloc(sizeof(unsigned int)
      *(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL));

  /* Some tings to try:
   *
   * 1) Terminate the 2D index with a leading "1" bit, like Warren/Salmon to
   * indicate the width of the key and therefore its tier.
   */

  /* Loop over A. */
  for (A_k_lookup_index = 0, B_k_lookup_index = 0; A_k_lookup_index < A_k_lookup.size-1; A_k_lookup_index++)
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
          number_dropped_blocks++;
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
          printf("[FIXME] C block missing\n");
          exit(1);
        }

        multiply_stream[stream_index].C = C_data;

        /* Done with this stream element. */
        stream_index++;
      }

      /* Test how quickly we tested out in the previous loop over B. */
      if (j == B_k_lookup.index[B_k_lookup_index])
      {
        /* We never went past the first block in B. Since k segments are norm
         * sorted, we know that we can skip the rest of this k segment in A. */
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

      /* Increment loop counter. */
      i++;
    }
  }

  /* Check. */
  if (stream_index > (A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL))
  {
    printf("multiply stream has too many elements, has %u but is only dimensioned for %u\n", stream_index,
        (A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL));
    exit(1);
  }

  spamm_timer_stop(timer);
  convolute_timer = spamm_timer_get(timer);
  printf("%llu timer units\n", convolute_timer);

  printf("[multiply] dropped %u blocks, placed %u blocks into stream\n", number_dropped_blocks, stream_index);
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
  free_timer = spamm_timer_get(timer);
  printf("%llu timer units\n", free_timer);
#endif

#ifdef SPAMM_MULTIPLY_STREAM
  /* Call stream product. */
  printf("[multiply] stream multiply ");
  spamm_timer_start(timer);

#if defined(SPAMM_KERNEL_IMPLEMENTATION_NO_CHECK)
  printf("(no-checks kernel)... ");
  spamm_stream_kernel_no_checks_SSE(stream_index, alpha, tolerance, multiply_stream);
#elif defined(SPAMM_KERNEL_IMPLEMENTATION_NO_CHECK_SSE4_1)
  printf("(no-checks kernel SSE4.1)... ");
  spamm_stream_kernel_no_checks_SSE4_1(stream_index, alpha, tolerance, multiply_stream);
#elif defined(SPAMM_KERNEL_IMPLEMENTATION_EXPERIMENTAL)
  printf("(experimental kernel)... ");
  spamm_stream_kernel_C(stream_index, alpha, tolerance, multiply_stream);
#elif defined(SPAMM_KERNEL_IMPLEMENTATION_STANDARD)
  printf("(standard kernel)... ");
  spamm_stream_kernel_SSE(stream_index, alpha, tolerance, multiply_stream);
#elif defined(SPAMM_KERNEL_IMPLEMENTATION_STANDARD_SSE4_1)
  printf("(standard kernel SSE4.1)... ");
  spamm_stream_kernel_SSE4_1(stream_index, alpha, tolerance, multiply_stream);
#elif defined(SPAMM_KERNEL_IMPLEMENTATION_Z_CURVE_SSE)
  printf("(Z-curve kernel SSE)... ");
  spamm_stream_kernel_Z_curve_SSE4_1(stream_index, alpha, tolerance, multiply_stream);
#elif defined(SPAMM_KERNEL_IMPLEMENTATION_Z_CURVE_SSE4_1)
  printf("(Z-curve kernel SSE4.1)... ");
  spamm_stream_kernel_Z_curve_SSE4_1(stream_index, alpha, tolerance, multiply_stream);
#else
  printf("(unknown kernel)... ");
  exit(1);
#endif

  spamm_timer_stop(timer);
  if (timer_type == papi_flop)
  {
    /* We multiply the floprate by 4 because most of the flops come from
     * single precision SSE instructions, each of which is counted only once
     * in the floprate. */
    printf("%1.2f Gflop/s, ", 4*spamm_timer_get_floprate(timer)/1e3);
  }
  stream_timer = spamm_timer_get(timer);
  printf("%llu timer units\n", stream_timer);
#endif

#ifdef SPAMM_MULTIPLY_PRODUCT_COUNT
  /* Loop through the stream and determine the number of products. */
  printf("[multiply] counting number of block products... ");
  fflush(stdout);
  for (index = 0; index < stream_index; index++)
  {
    for (i = 0; i < SPAMM_N_KERNEL_BLOCK; i++) {
      for (j = 0; j < SPAMM_N_KERNEL_BLOCK; j++) {
        for (k = 0; k < SPAMM_N_KERNEL_BLOCK; k++)
        {
          if (multiply_stream[index].A->norm[i*SPAMM_N_KERNEL_BLOCK+k]*multiply_stream[index].B->norm[k*SPAMM_N_KERNEL_BLOCK+j] >= tolerance)
          {
            number_products++;
          }
        }
      }
    }
  }
  printf("%u products out of %u possible (%1.2f%%)\n", number_products,
      (A->N/SPAMM_N_BLOCK)*(A->N/SPAMM_N_BLOCK)*(A->N/SPAMM_N_BLOCK),
      (double) number_products/(double) ((A->N/SPAMM_N_BLOCK)*(A->N/SPAMM_N_BLOCK)*(A->N/SPAMM_N_BLOCK))*100);
#endif

#ifdef SPAMM_MULTIPLY_FINAL_FREE
  /* Free memory. */
  free(multiply_stream);
  spamm_timer_delete(&timer);
#endif
}
