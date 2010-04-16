#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

double
spamm_get_element (const int i, const int j, const struct spamm_node_t *node)
{
  assert(node != NULL);

  int l, k;

  /* Test whether we are at a leaf node. */
  if (node->child == NULL)
  {
    if (node->block_dense != NULL)
    {
      return node->block_dense[spamm_dense_index(i-node->M_lower, j-node->N_lower, node->N_block)];
    }

    else { return 0.0; }
  }

  else
  {
    for (l = 0; l < node->M_child; ++l) {
      for (k = 0; k < node->N_child; ++k)
      {
        if (i >= (node->M_lower+(node->M_upper-node->M_lower)*l/node->M_child) &&
            i <  (node->M_lower+(node->M_upper-node->M_lower)*(l+1)/node->M_child) &&
            j >= (node->N_lower+(node->N_upper-node->N_lower)*k/node->N_child) &&
            j <  (node->N_lower+(node->N_upper-node->N_lower)*(k+1)/node->N_child))
        {
          return spamm_get_element(i, j, node->child[spamm_dense_index(l, k, node->N_child)]);
        }
      }
    }
  }
}

double
spamm_get (const int i, const int j, const struct spamm_t *A)
{
  assert(A != NULL);

  if (i < 0)     { spamm_log("(i = %i) < 0\n",  __FILE__, __LINE__, i); exit(1); }
  if (j < 0)     { spamm_log("(j = %i) < 0\n",  __FILE__, __LINE__, j); exit(1); }
  if (i >= A->M) { spamm_log("(i = %i) >= (M = %i)\n", __FILE__, __LINE__, i, A->M); exit(1); }
  if (j >= A->N) { spamm_log("(j = %i) >= (N = %i)\n", __FILE__, __LINE__, j, A->N); exit(1); }

  /* Recurse down to find the element. */
  if (A->root != NULL)
  {
    return spamm_get_element(i, j, A->root);
  }

  else { return 0.0; }
}
