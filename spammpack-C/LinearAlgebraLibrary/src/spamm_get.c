#include "spamm.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

double
spamm_get_element (const int i, const int j, const struct spamm_node_t *node)
{
  /* Test whether we are at a leaf node. */
  if (node->block_dense != NULL)
  {
    return node->block_dense[spamm_dense_index(i-node->child_M, j-node->child_N, node->block_N)];
  }

  else
  {
    if (i < node->child_M/2 && j < node->child_N/2) { return spamm_get_element(i, j, node->child[spamm_dense_index(0, 0, 2)]); }
    else if (i < node->child_M/2)                   { return spamm_get_element(i, j, node->child[spamm_dense_index(0, 1, 2)]); }
    else if (j < node->child_M/2)                   { return spamm_get_element(i, j, node->child[spamm_dense_index(1, 0, 2)]); }
    else                                            { return spamm_get_element(i, j, node->child[spamm_dense_index(1, 1, 2)]); }
  }
}

double
spamm_get (const int i, const int j, const struct spamm_t *A)
{
  assert(A != NULL);

  if (i < 0) { fprintf(stderr, "[%s:%i] i < 0\n", __FILE__, __LINE__); }
  if (j < 0) { fprintf(stderr, "[%s:%i] j < 0\n", __FILE__, __LINE__); }
  if (i >= A->M) { fprintf(stderr, "[%s:%i] i >= M\n", __FILE__, __LINE__); }
  if (i >= A->N) { fprintf(stderr, "[%s:%i] i >= N\n", __FILE__, __LINE__); }

  /* Recurse down to find the element. */
  if (A->root != NULL)
  {
    return spamm_get_element(i, j, A->root);
  }

  else { return 0.0; }
}
