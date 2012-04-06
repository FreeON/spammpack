#include "spamm.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Get an element from a matrix.
 *
 * @param i The row index.
 * @param j The column index.
 * @param A The matrix.
 *
 * @return The matrix element \f$A(i,j)\f$.
 */
float
spamm_get (const unsigned int i, const unsigned int j, const struct spamm_t *A)
{
  unsigned int index, i_tier, j_tier, delta_index;
  struct spamm_hashtable_t *node_hashtable;
  struct spamm_data_t *data;
  float Aij = 0;

  assert(A != NULL);

  if (i >= A->M || j >= A->N)
  {
    fprintf(stderr, "illegal index values for A_ij\n");
    exit(1);
  }

  /* Go into kernel tier hash and retrieve proper node. */
  delta_index = (unsigned int) floor(A->N_padded/pow(SPAMM_N_CHILD, A->kernel_tier));

  i_tier = i/delta_index;
  j_tier = j/delta_index;

  /* Construct linear index of the node on this tier. */
  index = spamm_index_2D(i_tier, j_tier);

  /* Get hash table at this tier. */
  node_hashtable = A->tier_hashtable[A->kernel_tier];

  if ((data = spamm_hashtable_lookup(node_hashtable, index)) != NULL)
  {
    Aij = data->block_dense[spamm_index_kernel_block(i%delta_index, j%delta_index, A->layout)];
  }

  return Aij;
}
