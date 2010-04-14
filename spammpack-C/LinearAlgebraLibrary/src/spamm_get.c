#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

double
spamm_get_element (const int i, const int j, const struct spamm_node_t *node)
{
  /* Test whether we are at a leaf node. */
  if (node->block_dense != NULL)
  {
    return node->block_dense[spamm_dense_index(i-node->M_lower, j-node->N_lower, node->N_block)];
  }

  else
  {
    if (i < (node->M_lower+(node->M_upper-node->M_lower)/2) && j < (node->N_lower+(node->N_upper-node->N_lower)/2))
    {
      return spamm_get_element(i, j, node->child[spamm_dense_index(0, 0, 2)]);
    }

    else if (i < (node->M_lower+(node->M_upper-node->M_lower)/2))
    {
      return spamm_get_element(i, j, node->child[spamm_dense_index(0, 1, 2)]);
    }

    else if (j < (node->N_lower+(node->N_upper-node->N_lower)/2))
    {
      return spamm_get_element(i, j, node->child[spamm_dense_index(1, 0, 2)]);
    }

    else
    {
      return spamm_get_element(i, j, node->child[spamm_dense_index(1, 1, 2)]);
    }
  }
}

double
spamm_get (const int i, const int j, const struct spamm_t *A)
{
  assert(A != NULL);

  if (i < 0)     { spamm_log("i < 0\n",  __FILE__, __LINE__); exit(1); }
  if (j < 0)     { spamm_log("j < 0\n",  __FILE__, __LINE__); exit(1); }
  if (i >= A->M) { spamm_log("i >= M\n", __FILE__, __LINE__); exit(1); }
  if (i >= A->N) { spamm_log("i >= N\n", __FILE__, __LINE__); exit(1); }

  /* Recurse down to find the element. */
  if (A->root != NULL)
  {
    return spamm_get_element(i, j, A->root);
  }

  else { return 0.0; }
}
