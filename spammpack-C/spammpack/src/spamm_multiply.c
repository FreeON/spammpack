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

  unsigned int x_masked = (*x) & 0x55555555;
  unsigned int y_masked = (*y) & 0x55555555;

  if (x_masked < y_masked)       { return -1; }
  else if (x_masked == y_masked) { return 0; }
  else                           { return 1; }
}

gint
spamm_multiply_compare_index_column (gconstpointer a, gconstpointer b)
{
  const unsigned int *x = a;
  const unsigned int *y = b;

  unsigned int x_masked = (*x) & 0xaaaaaaaa;
  unsigned int y_masked = (*y) & 0xaaaaaaaa;

  if (x_masked < y_masked)       { return -1; }
  else if (x_masked == y_masked) { return 0; }
  else                           { return 1; }
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
spamm_multiply (const float alpha, struct spamm_t *A, struct spamm_t *B,
    const float beta, struct spamm_t *C)
{
  GHashTable *A_tier_hashtable;
  GHashTable *B_tier_hashtable;
  GList *A_index_unsorted;
  GList *A_index_sorted;
  GList *B_index_unsorted;
  GList *B_index_sorted;
  struct spamm_multiply_index_list_t A_index;
  struct spamm_multiply_index_list_t B_index;

  struct spamm_data_t *data;

  unsigned int i, j;
  unsigned int first_B_index;
  unsigned int last_B_index;
  unsigned int last_A_index;
  unsigned int convolution_index;
  unsigned int k_match_index;

  char bitstring[1000];

  assert(A != NULL);
  assert(B != NULL);
  assert(C != NULL);

  /* Multiply C with beta. */

  /* Sort 2D indices on k, i.e. either on row or column index. */
  A_tier_hashtable = g_hash_table_lookup(A->tier_hashtable, &A->kernel_tier);
  A_index_unsorted = g_hash_table_get_keys(A_tier_hashtable);
  A_index_sorted = g_list_sort(A_index_unsorted, spamm_multiply_compare_index_column);

  B_tier_hashtable = g_hash_table_lookup(B->tier_hashtable, &B->kernel_tier);
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
  first_B_index = 0;
  last_A_index = (A_index.index_3D[0] & MASK_3D_K);
  last_B_index = B_index.last_index;
  for (i = 0; i < A_index.last_index; i++)
  {
    spamm_uint_to_bin_string(A_index.index_2D[i], bitstring);
    printf("A_index_2D[%u] = %u (%s), ", i, A_index.index_2D[i], bitstring);
    spamm_uint_to_bin_string(A_index.index_3D[i], bitstring);
    printf("A_index_3D[%u] = %u (%s), ", i, A_index.index_3D[i], bitstring);
    spamm_uint_to_bin_string(last_A_index, bitstring);
    printf("last_A_index = %u (%s)\n", last_A_index, bitstring);

    if (last_A_index != (A_index.index_3D[i] & MASK_3D_K))
    {
      printf("i = %u, new last_A_index %u --> ", i, last_A_index);
      last_A_index = (A_index.index_3D[i] & MASK_3D_K);
      first_B_index = last_B_index;
      last_B_index = B_index.last_index;
      printf("%u, B index starting from %u\n", last_A_index, first_B_index);
    }

    for (j = first_B_index; j < last_B_index; j++)
    {
      spamm_uint_to_bin_string(B_index.index_2D[j], bitstring);
      printf("B_index_2D[%u] = %u (%s), ", i, B_index.index_2D[j], bitstring);
      spamm_uint_to_bin_string(B_index.index_3D[j], bitstring);
      printf("B_index_3D[%u] = %u (%s), ", i, B_index.index_3D[j], bitstring);
      spamm_uint_to_bin_string(B_index.index_3D[j] & MASK_3D_K, bitstring);
      printf("k masked B_index_3D[%u] = %u (%s)\n", i, B_index.index_3D[j] & MASK_3D_K, bitstring);

      convolution_index = (A_index.index_3D[i] & MASK_3D_IJ) & (B_index.index_3D[j] & MASK_3D_IJ);
      k_match_index = last_A_index ^ (B_index.index_3D[j] & MASK_3D_K);

      if (k_match_index != 0)
      {
        /* The k indices do not match. */
        printf("last_B_index = %u\n", j);
        last_B_index = j;
        break;
      }

      printf("convolution: %x\n", convolution_index);
    }
  }

  /* Call stream product. */
}
