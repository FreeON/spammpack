#include "spamm.h"
#include <assert.h>

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
  unsigned int last_index;
  unsigned int *index_2D;
  unsigned int *index_3D;
};

gint
spamm_multiply_compare_index_row (gconstpointer a, gconstpointer b)
{
  const unsigned int *x = a;
  const unsigned int *y = b;

  unsigned int x_masked = (*x) & MASK_2D_I;
  unsigned int y_masked = (*y) & MASK_2D_I;

  if (x_masked < y_masked)       { return -1; }
  else if (x_masked == y_masked) { return  0; }
  else                           { return  1; }
}

gint
spamm_multiply_compare_index_column (gconstpointer a, gconstpointer b)
{
  const unsigned int *x = a;
  const unsigned int *y = b;

  unsigned int x_masked = (*x) & MASK_2D_J;
  unsigned int y_masked = (*y) & MASK_2D_J;

  if (x_masked < y_masked)       { return -1; }
  else if (x_masked == y_masked) { return  0; }
  else                           { return  1; }
}

void
spamm_multiply_copy_to_array (gpointer data, gpointer user_data)
{
  unsigned int *index = data;
  struct spamm_multiply_index_list_t *A_index = user_data;

  A_index->index_2D[A_index->last_index] = *index;
  (A_index->last_index)++;
}

void
spamm_multiply_beta (const float beta, struct spamm_t *A)
{
}

void
spamm_multiply (const float tolerance,
    const float alpha, struct spamm_t *A, struct spamm_t *B,
    const float beta, struct spamm_t *C)
{
  GHashTable *A_tier_hashtable;
  GHashTable *B_tier_hashtable;
  GHashTable *C_tier_hashtable;
  GList *A_index_unsorted;
  GList *A_index_sorted;
  GList *B_index_unsorted;
  GList *B_index_sorted;
  struct spamm_multiply_index_list_t A_index;
  struct spamm_multiply_index_list_t B_index;

  struct spamm_data_t *A_block;
  struct spamm_data_t *B_block;
  struct spamm_data_t *C_block;

  struct spamm_data_t *data;

  unsigned int i, j, k;
  unsigned int first_B_index;
  unsigned int last_B_index;
  unsigned int last_A_index;
  unsigned int convolution_index;
  unsigned int convolution_index_2D;
  unsigned int k_match_index;
  unsigned int C_i, C_j;

  struct multiply_stream_t *multiply_stream;
  unsigned int stream_index;

  char bitstring[1000];

  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);

  /* Multiply C with beta. */
  spamm_multiply_beta(beta, C);

  /* Sort 2D indices on k, i.e. either on row or column index. */
  A_tier_hashtable = g_hash_table_lookup(A->tier_hashtable, &A->kernel_tier);
  B_tier_hashtable = g_hash_table_lookup(B->tier_hashtable, &B->kernel_tier);
  C_tier_hashtable = g_hash_table_lookup(C->tier_hashtable, &C->kernel_tier);

  A_index_unsorted = g_hash_table_get_keys(A_tier_hashtable);
  A_index_sorted = g_list_sort(A_index_unsorted, spamm_multiply_compare_index_column);

  B_index_unsorted = g_hash_table_get_keys(B_tier_hashtable);
  B_index_sorted = g_list_sort(B_index_unsorted, spamm_multiply_compare_index_row);

  /* Copy sorted indices to array for quick access. */
  A_index.last_index = 0;
  A_index.index_2D = (unsigned int*) malloc(sizeof(unsigned int)*g_list_length(A_index_sorted));
  A_index.index_3D = (unsigned int*) malloc(sizeof(unsigned int)*g_list_length(A_index_sorted));

  B_index.last_index = 0;
  B_index.index_2D = (unsigned int*) malloc(sizeof(unsigned int)*g_list_length(B_index_sorted));
  B_index.index_3D = (unsigned int*) malloc(sizeof(unsigned int)*g_list_length(B_index_sorted));

  g_list_foreach(A_index_sorted, spamm_multiply_copy_to_array, &A_index);
  g_list_foreach(B_index_sorted, spamm_multiply_copy_to_array, &B_index);

  /* Copy appropriate 3D convolution index to arrays. */
  for (i = 0; i < A_index.last_index; i++)
  {
    data = g_hash_table_lookup(A_tier_hashtable, &A_index.index_2D[i]);
    A_index.index_3D[i] = data->index_3D_ik0;
  }

  for (i = 0; i < B_index.last_index; i++)
  {
    data = g_hash_table_lookup(B_tier_hashtable, &B_index.index_2D[i]);
    B_index.index_3D[i] = data->index_3D_0kj;
  }

  /* Free some memory. */
  g_list_free(A_index_unsorted);
  g_list_free(B_index_unsorted);

  g_list_free(A_index_sorted);
  g_list_free(B_index_sorted);

  /* Convolute by constructing product 3D index. */
  stream_index = 0;
  multiply_stream = (struct multiply_stream_t*) malloc(sizeof(struct multiply_stream_t)*(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL)*(A->N_padded/SPAMM_N_KERNEL));

  first_B_index = 0;
  last_A_index = (A_index.index_3D[0] & MASK_3D_K);
  last_B_index = B_index.last_index;
  for (i = 0; i < A_index.last_index; i++)
  {
    if (last_A_index != (A_index.index_3D[i] & MASK_3D_K))
    {
      last_A_index = (A_index.index_3D[i] & MASK_3D_K);
      first_B_index = last_B_index;
      last_B_index = B_index.last_index;
    }

    /* Get reference to dense block of A. */
    A_block = g_hash_table_lookup(A_tier_hashtable, &A_index.index_2D[i]);

    for (j = first_B_index; j < last_B_index; j++)
    {
      /* Get reference to dense block of B. */
      B_block = g_hash_table_lookup(B_tier_hashtable, &B_index.index_2D[j]);

      convolution_index = (A_index.index_3D[i] & MASK_3D_IJ) | (B_index.index_3D[j] & MASK_3D_IJ);
      k_match_index = last_A_index ^ (B_index.index_3D[j] & MASK_3D_K);

      if (k_match_index != 0)
      {
        /* The k indices do not match. */
        last_B_index = j;
        break;
      }

      convolution_index_2D = spamm_index_3D_i0j_to_2D(convolution_index);

      /* Get reference to dense block of C. */
      C_block = g_hash_table_lookup(C_tier_hashtable, &convolution_index_2D);

      multiply_stream[stream_index].A_block = A_block->block_dense_dilated;
      multiply_stream[stream_index].B_block = B_block->block_dense;
      multiply_stream[stream_index].C_block = C_block->block_dense;

      for (k = 0; k < 16; k++)
      {
        multiply_stream[stream_index].norm[k] = A_block->norm[k];
      }

      for (k = 16; k < 32; k++)
      {
        multiply_stream[stream_index].norm[k] = B_block->norm[k-16];
      }

      stream_index++;
    }
  }
  free(A_index.index_2D);
  free(A_index.index_3D);
  free(B_index.index_2D);
  free(B_index.index_3D);

  /* Call stream product. */
  spamm_stream_kernel(stream_index, alpha, tolerance, multiply_stream);

  /* Free memory. */
  free(multiply_stream);
}
