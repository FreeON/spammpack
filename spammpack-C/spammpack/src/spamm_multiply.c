#include "spamm.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

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
#define MASK_3D_K  0x12492492
#define MASK_3D_IJ 0x2DB6DB6D
#define MASK_2D_J  0x55555555
#define MASK_2D_I  0xaaaaaaaa

struct spamm_multiply_index_list_t
{
  unsigned int size;
  unsigned int *index_2D;
  unsigned int *index_3D;
};

struct spamm_multiply_k_lookup_t
{
  unsigned int size;
  unsigned int *index;
};

gint
spamm_multiply_compare_index_row (gconstpointer a, gconstpointer b, gpointer user_data)
{
  const unsigned int *x = a;
  const unsigned int *y = b;

  GHashTable *tier_hashtable = user_data;

  struct spamm_data_t *a_data;
  struct spamm_data_t *b_data;

  unsigned int x_masked = (*x) & MASK_2D_I;
  unsigned int y_masked = (*y) & MASK_2D_I;

  if (x_masked < y_masked)      { return -1; }
  else if (x_masked > y_masked) { return  1; }
  else
  {
    /* Compare norms. */
    a_data = g_hash_table_lookup(tier_hashtable, x);
    b_data = g_hash_table_lookup(tier_hashtable, y);

    /* Sort norms within k in descending order. */
    if (a_data->norm > b_data->norm)      { return -1; }
    else if (a_data->norm < b_data->norm) { return  1; }
    else                                  { return  0; }
  }
}

gint
spamm_multiply_compare_index_column (gconstpointer a, gconstpointer b, gpointer user_data)
{
  const unsigned int *x = a;
  const unsigned int *y = b;

  GHashTable *tier_hashtable = user_data;

  struct spamm_data_t *a_data;
  struct spamm_data_t *b_data;

  unsigned int x_masked = (*x) & MASK_2D_J;
  unsigned int y_masked = (*y) & MASK_2D_J;

  if (x_masked < y_masked)       { return -1; }
  else if (x_masked > y_masked)  { return  1; }
  else
  {
    /* Compare norms. */
    a_data = g_hash_table_lookup(tier_hashtable, x);
    b_data = g_hash_table_lookup(tier_hashtable, y);

    /* Sort norms within k in descending order. */
    if (a_data->norm > b_data->norm)      { return -1; }
    else if (a_data->norm < b_data->norm) { return  1; }
    else                                  { return  0; }
  }
}

void
spamm_multiply_copy_to_array (gpointer data, gpointer user_data)
{
  unsigned int *index = data;
  struct spamm_multiply_index_list_t *A_index = user_data;

  A_index->index_2D[A_index->size] = *index;
  (A_index->size)++;
}

void
spamm_multiply_beta_block (gpointer key, gpointer value, gpointer user_data)
{
  struct spamm_data_t *data = value;
  float *beta = user_data;
  unsigned int i;

  for (i = 0; i < SPAMM_N_KERNEL*SPAMM_N_KERNEL; i++)
  {
    data->block_dense[i] *= (*beta);

    data->block_dense_dilated[i+0] *= (*beta);
    data->block_dense_dilated[i+1] *= (*beta);
    data->block_dense_dilated[i+2] *= (*beta);
    data->block_dense_dilated[i+3] *= (*beta);
  }
}

void
spamm_multiply_beta (const float beta, struct spamm_t *A)
{
  GHashTable *tier_hashtable;

  tier_hashtable = g_hash_table_lookup(A->tier_hashtable, &A->kernel_tier);
  g_hash_table_foreach(tier_hashtable, spamm_multiply_beta_block, (void*) &beta);
}

/** Multiply to matrices, i.e. \f$ C = \alpha A \times B + \beta C\f$.
 *
 * @param tolerance The SpAMM tolerance of this product.
 * @param alpha The paramater \f$\alpha\f$.
 * @param A The matrix \f$A\f$.
 * @param B The matrix \f$B\f$.
 * @param beta The paramater \f$\beta\f$.
 * @param C The matrix \f$C\f$.
 */
void
spamm_multiply (const float tolerance,
    const float alpha, struct spamm_t *A, struct spamm_t *B,
    const float beta, struct spamm_t *C)
{
  GHashTable *A_tier_hashtable;
  GHashTable *B_tier_hashtable;
  GHashTable *C_tier_hashtable;

  GList *A_index_sorted;
  GList *B_index_sorted;
  GList *list_element;

  struct spamm_multiply_index_list_t A_index;
  struct spamm_multiply_index_list_t B_index;

  struct spamm_data_t *A_block;
  struct spamm_data_t *B_block;
  struct spamm_data_t *C_block;

  struct spamm_data_t *data;

  unsigned int i, j, k, k_check;
  unsigned int index;
  unsigned int convolution_index;
  unsigned int convolution_index_2D;
  unsigned int A_k_lookup_index;
  unsigned int B_k_lookup_index;
  unsigned int A_k, B_k;

  struct spamm_multiply_k_lookup_t A_k_lookup;
  struct spamm_multiply_k_lookup_t B_k_lookup;

  struct multiply_stream_t *multiply_stream;
  unsigned int stream_index;

  struct spamm_timer_t *beta_timer          = spamm_timer_new();
  struct spamm_timer_t *sort_timer          = spamm_timer_new();
  struct spamm_timer_t *k_lookuptable_timer = spamm_timer_new();
  struct spamm_timer_t *copy_timer          = spamm_timer_new();
  struct spamm_timer_t *copy_3D_timer       = spamm_timer_new();
  struct spamm_timer_t *free_timer          = spamm_timer_new();
  struct spamm_timer_t *convolute_timer     = spamm_timer_new();
  struct spamm_timer_t *stream_timer        = spamm_timer_new();
  struct spamm_timer_t *free_2_timer        = spamm_timer_new();

  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);

  /* Multiply C with beta. */
  printf("[multiply] multiplying C with beta... ");
  spamm_timer_start(beta_timer);

  spamm_multiply_beta(beta, C);

  spamm_timer_stop(beta_timer);
  printf("%1.2e s\n", spamm_timer_get_seconds(beta_timer));

  /* Sort 2D indices on k, i.e. either on row or column index. */
  printf("[multiply] sorting A and B... ");
  spamm_timer_start(sort_timer);

  A_tier_hashtable = g_hash_table_lookup(A->tier_hashtable, &A->kernel_tier);
  B_tier_hashtable = g_hash_table_lookup(B->tier_hashtable, &B->kernel_tier);
  C_tier_hashtable = g_hash_table_lookup(C->tier_hashtable, &C->kernel_tier);

  A_index_sorted = g_hash_table_get_keys(A_tier_hashtable);
  A_index_sorted = g_list_sort_with_data(A_index_sorted, spamm_multiply_compare_index_column, A_tier_hashtable);

  B_index_sorted = g_hash_table_get_keys(B_tier_hashtable);
  B_index_sorted = g_list_sort_with_data(B_index_sorted, spamm_multiply_compare_index_row, B_tier_hashtable);

  spamm_timer_stop(sort_timer);
  printf("%1.2e s\n", spamm_timer_get_seconds(sort_timer));

  /* Create a lookup table for the start of a particular k index in the sorted
   * arrays. */
  printf("[multiply] creating k lookup trables... ");
  spamm_timer_start(k_lookuptable_timer);

  A_k_lookup.index = (unsigned int*) malloc(sizeof(unsigned int)*(A->N_padded/SPAMM_N_KERNEL+1));
  B_k_lookup.index = (unsigned int*) malloc(sizeof(unsigned int)*(B->N_padded/SPAMM_N_KERNEL+1));

  /* The index in A_k_lookup. */
  A_k_lookup.size = 0;

  /* The index in A_index_sorted. */
  i = 0;

  /* The last k value. Initially place it behind the largest expected k value. */
  k = A->N+1;

  for (list_element = g_list_first(A_index_sorted); list_element != NULL; list_element = g_list_next(list_element))
  {
    /* Extract k index from linear index. */
    index = *((unsigned int*) list_element->data);
    spamm_index_2D_to_ij(index, NULL, &k_check);

    if (k != k_check)
    {
      A_k_lookup.index[A_k_lookup.size++] = i;
      k = k_check;
    }

    i++;
  }

  /* Add terminating entry to lookup list. */
  A_k_lookup.index[A_k_lookup.size++] = g_list_length(A_index_sorted);

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

  for (list_element = g_list_first(B_index_sorted); list_element != NULL; list_element = g_list_next(list_element))
  {
    /* Extract k index from linear index. */
    index = *((unsigned int*) list_element->data);
    spamm_index_2D_to_ij(index, &k_check, NULL);

    if (k != k_check)
    {
      B_k_lookup.index[B_k_lookup.size++] = i;
      k = k_check;
    }

    i++;
  }

  /* Add terminating entry to lookup list. */
  B_k_lookup.index[B_k_lookup.size++] = g_list_length(B_index_sorted);

  /* Check. */
  if (B_k_lookup.size > B->N_padded/SPAMM_N_KERNEL+1)
  {
    printf("k lookup table too long for B, estimated %u elemens, but found %u\n", B->N_padded/SPAMM_N_KERNEL+1, B_k_lookup.size);
    exit(1);
  }

  spamm_timer_stop(k_lookuptable_timer);
  printf("%1.2e s\n", spamm_timer_get_seconds(k_lookuptable_timer));

  /* Copy sorted indices to array for quick access. */
  printf("[multiply] copying indices to array... ");
  spamm_timer_start(copy_timer);

  A_index.size = 0;
  A_index.index_2D = (unsigned int*) malloc(sizeof(unsigned int)*g_list_length(A_index_sorted));
  A_index.index_3D = (unsigned int*) malloc(sizeof(unsigned int)*g_list_length(A_index_sorted));

  B_index.size = 0;
  B_index.index_2D = (unsigned int*) malloc(sizeof(unsigned int)*g_list_length(B_index_sorted));
  B_index.index_3D = (unsigned int*) malloc(sizeof(unsigned int)*g_list_length(B_index_sorted));

  g_list_foreach(A_index_sorted, spamm_multiply_copy_to_array, &A_index);
  g_list_foreach(B_index_sorted, spamm_multiply_copy_to_array, &B_index);

  spamm_timer_stop(copy_timer);
  printf("%1.2e s\n", spamm_timer_get_seconds(copy_timer));

  /* Copy appropriate 3D convolution index to arrays. */
  printf("[multiply] copying 3D convolution index to arrays... ");
  spamm_timer_start(copy_3D_timer);

  for (i = 0; i < A_index.size; i++)
  {
    data = g_hash_table_lookup(A_tier_hashtable, &A_index.index_2D[i]);
    A_index.index_3D[i] = data->index_3D_ik0;
  }

  for (i = 0; i < B_index.size; i++)
  {
    data = g_hash_table_lookup(B_tier_hashtable, &B_index.index_2D[i]);
    B_index.index_3D[i] = data->index_3D_0kj;
  }

  spamm_timer_stop(copy_3D_timer);
  printf("%1.2e s\n", spamm_timer_get_seconds(copy_3D_timer));

  /* Free some memory. */
  printf("[multiply] free some memory... ");
  spamm_timer_start(free_timer);

  g_list_free(A_index_sorted);
  g_list_free(B_index_sorted);

  spamm_timer_stop(free_timer);
  printf("%1.2e s\n", spamm_timer_get_seconds(free_timer));

  /* Convolute by constructing product 3D index. */
  printf("[multiply] convolute... ");
  spamm_timer_start(convolute_timer);

  stream_index = 0;
  multiply_stream = (struct multiply_stream_t*) malloc(sizeof(struct multiply_stream_t)*(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL));

  /* Loop over A. */
  A_k_lookup_index = 0;
  B_k_lookup_index = 0;
  for (i = 0; i < A_index.size; )
  {
    /* Get k value of A. */
    A_k = spamm_index_3D_ikj_to_k(A_index.index_3D[i]);

    /* Get k value of B. */
    B_k = spamm_index_3D_ikj_to_k(B_index.index_3D[B_k_lookup.index[B_k_lookup_index]]);

    /* Compare k values. */
    if (A_k > B_k)
    {
      /* Advance B in k. */
      B_k_lookup_index++;
      continue;
    }

    else if (A_k < B_k)
    {
      /* Advance A in k. */
      A_k_lookup_index++;

      /* Set loop counter correctly. */
      i = A_k_lookup.index[A_k_lookup_index];
      continue;
    }

    /* Get reference to dense block of A. */
    A_block = g_hash_table_lookup(A_tier_hashtable, &A_index.index_2D[i]);

    /* Loop over subset of B with matching k. */
    for (j = B_k_lookup.index[B_k_lookup_index]; j < B_k_lookup.index[B_k_lookup_index+1]; j++)
    {
      /* Get reference to dense block of B. */
      B_block = g_hash_table_lookup(B_tier_hashtable, &B_index.index_2D[j]);

      /* Perform norm product and test whether to keep this term. */
      if (A_block->node_norm*B_block->node_norm <= tolerance)
      {
        break;
      }

      /* Get the linear 2D index of the C block. */
      convolution_index = (A_index.index_3D[i] & MASK_3D_IJ) | (B_index.index_3D[j] & MASK_3D_IJ);
      convolution_index_2D = spamm_index_3D_i0j_to_2D(convolution_index);

      /* Get reference to dense block of C. */
      C_block = g_hash_table_lookup(C_tier_hashtable, &convolution_index_2D);

      /* Set references to matrix block in multiply stream. */
      multiply_stream[stream_index].A_block = A_block->block_dense_dilated;
      multiply_stream[stream_index].B_block = B_block->block_dense;
      multiply_stream[stream_index].C_block = C_block->block_dense;

      /* Set the kernel block norms. */
      for (k = 0; k < 16; k++)
      {
        multiply_stream[stream_index].norm[k] = A_block->norm[k];
      }

      for (k = 16; k < 32; k++)
      {
        multiply_stream[stream_index].norm[k] = B_block->norm[k-16];
      }

      /* Done with this stream element. */
      stream_index++;
    }

    /* Test how quickly we tested out in the previous loop over B. */
    if (j == B_k_lookup.index[B_k_lookup_index])
    {
      /* We never went past the first block in B. Since k segments are norm
       * sorted, we know that we can skip the rest of this k segment in A. */
      A_k_lookup_index++;

      /* Set loop counter correctly. */
      i = A_k_lookup.index[A_k_lookup_index];
      continue;
    }

    /* Increment loop counter. */
    i++;
  }

  /* Check. */
  if (stream_index > (A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL))
  {
    printf("multiply stream has too many elements, has %u but is only dimensioned for %u\n", stream_index,
        (A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL));
    exit(1);
  }

  spamm_timer_stop(convolute_timer);
  printf("%1.2e s\n", spamm_timer_get_seconds(convolute_timer));

  /* Free memory. */
  printf("[multiply] free some more memory... ");
  spamm_timer_start(free_2_timer);

  free(A_index.index_2D);
  free(A_index.index_3D);
  free(B_index.index_2D);
  free(B_index.index_3D);

  spamm_timer_stop(free_2_timer);
  printf("%1.2e s\n", spamm_timer_get_seconds(free_2_timer));

  /* Call stream product. */
  printf("[multiply] stream multiply... ");
  spamm_timer_start(stream_timer);

  spamm_stream_kernel(stream_index, alpha, tolerance, multiply_stream);

  spamm_timer_stop(stream_timer);
  printf("%1.2e s\n", spamm_timer_get_seconds(stream_timer));

  /* Print out total time. */
  printf("[multiply] total time elapsed for everything but stream: %1.2e s\n",
      spamm_timer_get_seconds(beta_timer) +
      spamm_timer_get_seconds(sort_timer) +
      spamm_timer_get_seconds(copy_timer) +
      spamm_timer_get_seconds(copy_3D_timer) +
      spamm_timer_get_seconds(free_timer) +
      spamm_timer_get_seconds(convolute_timer));

  printf("[multiply] total time elapsed: %1.2e s\n",
      spamm_timer_get_seconds(beta_timer) +
      spamm_timer_get_seconds(sort_timer) +
      spamm_timer_get_seconds(copy_timer) +
      spamm_timer_get_seconds(copy_3D_timer) +
      spamm_timer_get_seconds(free_timer) +
      spamm_timer_get_seconds(convolute_timer) +
      spamm_timer_get_seconds(stream_timer));

  /* Free memory. */
  free(multiply_stream);

  spamm_timer_delete(&beta_timer);
  spamm_timer_delete(&sort_timer);
  spamm_timer_delete(&copy_timer);
  spamm_timer_delete(&copy_3D_timer);
  spamm_timer_delete(&free_timer);
  spamm_timer_delete(&convolute_timer);
  spamm_timer_delete(&stream_timer);
}
