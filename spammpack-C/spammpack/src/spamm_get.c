#include "spamm.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

float
spamm_get (const unsigned int i, const unsigned int j, const struct spamm_t *A)
{
  unsigned int index, i_tier, j_tier, delta_index;
  GHashTable *node_hashtable;
  struct spamm_data_t *data;
  float Aij = 0;

  assert(A != NULL);

  if (i >= A->M || j >= A->N)
  {
    printf("illegal index values for A_ij\n");
    exit(1);
  }

  /* Go into kernel tier hash and retrieve proper node. */
  delta_index = (unsigned int) floor(A->N_padded/pow(SPAMM_N_CHILD, A->kernel_tier));

  i_tier = (unsigned int) floor(i/delta_index);
  j_tier = (unsigned int) floor(j/delta_index);

  /* Construct linear index of the node on this tier. */
  index = spamm_index_2D(i_tier, j_tier);

  /* Get hash table at this tier. */
  node_hashtable = g_hash_table_lookup(A->tier_hashtable, &A->kernel_tier);

  if ((data = g_hash_table_lookup(node_hashtable, &index)) != NULL)
  {
    Aij = data->block_dense[i-i_tier*delta_index][j-j_tier*delta_index];
  }

  return Aij;
}
